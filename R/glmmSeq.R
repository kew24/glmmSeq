setClassUnion("character_or_list", c("character", "list"))
setClassUnion("df_or_matrix", c("data.frame", "matrix"))

#' An S4 class to define the glmmSeq output
#'
#' @slot formula The model formula
#' @slot stats the statistics from the glmm fit
#' @slot predict The predicted interception values
#' @slot reducedFormula The reduced formula with removed random effects
#' @slot countdata The input expression data
#' @slot metadata The input metadata
#' @slot modelData the model data for the glmm
#' @slot optInfo Information on whether the model was singular or converged
#' @slot errors Any errors
#' @slot variables The variables used in the formula

setClass("GlmmSeq", slots = list(
  formula = "formula",
  stats = "df_or_matrix",
  predict = "df_or_matrix",
  reducedFormula = "formula",
  countdata = "df_or_matrix",
  metadata = "df_or_matrix",
  modelData = "df_or_matrix",
  optInfo = "matrix",
  errors = "character_or_list",
  variables = "character_or_list"
))


#' Glmm for sequencing results
#'
#' @param modelFormula the model formula. For more information of formula
#' structure see \code{\link[lme4:glmer]{glmer}}
#' @param countdata the sequencing count data
#' @param metadata a data frame of sample information
#' @param id Column name in metadata which contains the sample IDs to be used
#' in pairing samples
#' @param dispersion a numeric vector of gene dispersion for the family. If
#' family does not use dispersion, this can be set as NA. 
#' @param sizeFactors size factors (default = NULL). If provided the glmer 
#' offset is set to log(sizeFactors). For more information see
#'  \code{\link[lme4:glmer]{glmer}}
#' @param reducedFormula Reduced design formula (default = "")
#' @param modelData Expanded design matrix
#' @param glmOnly Logical whether to only run 
#' \code{\link[stats:glm]{stats::glm}} (TRUE) instead of 
#' \code{\link[lme4:glmer]{lme4::glmer}} (FALSE). Default is FALSE. 
#' @param control the glmer control (default = glmerControl(optimizer = 
#' "bobyqa")). For more information see
#' \code{\link[lme4:glmerControl]{lme4::glmerControl}}.
#' @param family The GLM family, see 
#' \code{\link[stats:glm]{stats::glm}} and 
#' \code{\link[stats:family]{stats::family}}. If NULL 
#' \code{\link[MASS:negative.binomial]{MASS::negative.binomial}} is used. 
#' @param cores number of cores to use. Default = 1. 
#' @param removeDuplicatedMeasures whether to remove duplicated
#' conditions/repeated measurements for a given time point (default = FALSE).
#' @param removeSingles whether to remove individuals with only one measurement
#' (default = FALSE)
#' @param zeroCount numerical value to offset zeroes for the purpose of log
#' (default = 0.125)
#' @param verbose Logical whether to display messaging (default = TRUE)
#' @param returnList Logical whether to return results as a list or glmmSeq 
#' object (default = FALSE).
#' @param progress Logical whether to display a progress bar
#' @param ... Other parameters to pass to
#' \code{\link[lme4:glmer]{lme4::glmer()}}
#' @return Returns a GlmmSeq object with results for gene-wise general linear
#' mixed models or a list of results if returnList is TRUE.
#' @importFrom MASS negative.binomial
#' @importFrom lme4 subbars findbars glmer fixef glmerControl nobars isSingular
#' @importFrom stats update.formula model.matrix predict setNames
#' @importFrom parallel mclapply detectCores parLapply makeCluster clusterEvalQ
#' clusterExport stopCluster
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#' @importFrom car Anova
#' @importFrom methods slot new
#' @importFrom stats AIC complete.cases logLik reshape terms vcov
#' @export
#' @examples
#' data(PEAC_minimal_load)
#' disp <- apply(tpm, 1, function(x) {
#' (var(x, na.rm = TRUE)-mean(x, na.rm = TRUE))/(mean(x, na.rm = TRUE)**2)
#' })
#' MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      id = "PATID",
#'                      countdata = tpm["MS4A1", ],
#'                      metadata = metadata,
#'                      dispersion = disp["MS4A1"],
#'                      verbose = FALSE)
#' names(attributes(MS4A1glmm))


glmmSeq <- function(modelFormula,
                    countdata,
                    metadata,
                    id,
                    dispersion,
                    sizeFactors = NULL,
                    reducedFormula = "",
                    modelData = NULL,
                    glmOnly = FALSE,
                    control = glmerControl(optimizer = "bobyqa"),
                    family = NULL, 
                    cores = 1,
                    removeDuplicatedMeasures = FALSE,
                    removeSingles = FALSE,
                    zeroCount = 0.125,
                    verbose = TRUE,
                    returnList = FALSE, 
                    progress = TRUE,
                    ...) {
  
  # Catch errors
  if (! glmOnly & length(findbars(modelFormula)) == 0) {
    warning(paste("glmOnly is FALSE but no random effects specified in", 
                  "formula. Running glm instead"))
    glmOnly <- TRUE
  } 
  if (glmOnly & length(findbars(modelFormula)) > 0) {
    warning(paste("glmOnly is TRUE but random effects specified in formula.", 
                  "Running glmer instead"))
    glmOnly <- FALSE
  }
  if (ncol(countdata) != nrow(metadata)) {
    stop("countdata columns different size to metadata rows")
  }
  if (!is.null(sizeFactors) & ncol(countdata) != length(sizeFactors)) {
    stop("Different sizeFactors length")
  }
  if(is.null(family)){
    if (! all(rownames(countdata) %in% names(dispersion), nrow(countdata))) {
      stop("Dispersion length must match nrow in countdata")
    }
    family <- MASS::negative.binomial(theta = 1/dispersion)
  }
  if (! is.numeric(zeroCount)) stop("zeroCount must be numeric")
  if (zeroCount < 0) stop("zeroCount must be > = 0")
  if (zeroCount > 0) countdata[countdata == 0] <- zeroCount
  
  # Manipulate formulae
  fullFormula <- update.formula(modelFormula, count ~ ., simplify = FALSE)
  nonRandomFormula <- subbars(modelFormula)
  variables <- rownames(attr(terms(nonRandomFormula), "factors"))
  subsetMetadata <- metadata[, variables]
  ids <- as.character(metadata[, id])
  
  
  # Option to subset to remove duplicated timepoints
  if (removeDuplicatedMeasures) {
    # Check the distribution for duplicates
    check <- data.frame(table(droplevels(subsetMetadata)))
    check <- check[! check$Freq %in% c(0, 1), ]
    if (nrow(check) > 0) {
      mCheck <- as.character(apply(subsetMetadata[, variables], 1, function(x) {
        paste(as.character(x), collapse = " ")
      }))
      cCheck <- as.character(apply(check[, variables], 1, function(x) {
        paste(as.character(x), collapse = " ")
      }))
      countdata <- countdata[, ! mCheck %in% cCheck]
      sizeFactors <- sizeFactors[! mCheck %in% cCheck]
      subsetMetadata <- subsetMetadata[! mCheck %in% cCheck, ]
      ids <- droplevels(subsetMetadata[, id])
      warning(paste0(paste(check[, id], collapse = ", "),
                     " has multiple entries for identical ",
                     paste0(colnames(check)[! colnames(check) %in%
                                              c(id, "Freq")],
                            collapse = " and "),
                     ". These will all be removed."))
    }
  }
  
  
  # Option to subset to remove unpaired samples
  if (removeSingles) {
    singles <- names(table(ids)[table(ids) %in% c(0, 1)])
    nonSingleIDs <- which(! subsetMetadata[, id] %in% singles)
    
    countdata <- countdata[, nonSingleIDs]
    sizeFactors <- sizeFactors[nonSingleIDs]
    subsetMetadata <- subsetMetadata[nonSingleIDs, ]
    ids <- droplevels(subsetMetadata[, id])
  }
  
  # Check numbers and alignment
  if (! all(vapply(list(length(ids), nrow(subsetMetadata)), FUN = identical,
                   FUN.VALUE = TRUE, ncol(countdata)))) {
    stop("Alignment error: metadata rownames must match countdata colnames")
  }
  
  if (!is.null(sizeFactors)) offset <- log(sizeFactors) else offset <- NULL
  if (verbose) cat(paste0("\nn = ", length(ids), " samples, ",
                          length(unique(ids)), " individuals\n"))
  
  
  # setup model prediction
  if (reducedFormula == "") reducedFormula <- nobars(modelFormula)
  if (is.null(modelData)) {
    reducedVars <- rownames(attr(terms(reducedFormula), "factors"))
    varLevels <- lapply(reducedVars, function(x) {
      if (class(metadata[, x]) == "factor") {
        return(levels(subsetMetadata[, x]))
      } else {sort(unique(subsetMetadata[, x]))}
    })
    modelData <- expand.grid(varLevels)
    colnames(modelData) <- reducedVars
  }
  designMatrix <- model.matrix(reducedFormula, modelData)
  
  start <- Sys.time()
  fullList <- lapply(rownames(countdata), function(i) {
    familyInput <- family 
    if(length(familyInput$family) > 1) {
      familyInput$family <- familyInput$family[
        which(rownames(countdata) == i)]
    }
    list(y = countdata[i, ], family = familyInput)
  })
  
  # For each gene perform a fit
  if(glmOnly) {
    functionUsed <- glmApply 
    if("glmerControl" %in% class(control)) {
      control <- list()
      message("glmer control passed into glm, using control <- list() instead")
    }
  } else {
    functionUsed <- glmerApply
    if(! "glmerControl" %in% class(control)) {
      control <- glmerControl()
      message(paste("Control not of class 'glmer control', using control <-", 
                    "glmerControl() instead"))
    }
  }
  
  if (Sys.info()["sysname"] == "Windows") {
    cl <- makeCluster(cores)
    clusterExport(cl, varlist = c("functionUsed", 
                                  "fullList", "fullFormula",
                                  "subsetMetadata", "control", "modelData",
                                  "offset", "designMatrix", ...),
                  envir = environment())
    if (progress) {
      resultList <- pblapply(fullList, function(geneList) {
        functionUsed(geneList, fullFormula = fullFormula, data = subsetMetadata,
                     control = control, modelData = modelData, offset = offset,
                     designMatrix = designMatrix, ...)
      }, cl = cl)
    } else {
      resultList <- parLapply(cl = cl, fullList, function(geneList) {
        functionUsed(geneList, fullFormula = fullFormula, data = subsetMetadata,
                     control = control, modelData = modelData, offset = offset,
                     designMatrix = designMatrix, ...)
      })
    }
    stopCluster(cl)
  } else{
    if (progress) {
      resultList <- pbmclapply(fullList, function(geneList) {
        functionUsed(geneList, fullFormula = fullFormula, data = subsetMetadata,
                     control = control, modelData = modelData, offset = offset,
                     designMatrix = designMatrix)#, ...)
      }, mc.cores = cores)
    } else {
      resultList <- mclapply(fullList, function(geneList) {
        functionUsed(geneList, fullFormula = fullFormula, data = subsetMetadata,
                     control = control, modelData = modelData, offset = offset,
                     designMatrix = designMatrix, ...)
      }, mc.cores = cores)
    }
  }
  if ("value" %in% names(resultList)) {
    warning(paste0("Some warnings arising in generalized linear models: '", 
                   trimws(resultList$warning), "'"))
    resultList <- resultList$value
  }
  if(returnList) return(resultList)
  
  # Print timing if verbose
  time <- Sys.time() - start
  if (verbose) cat("\nCompleted in", time, units(time), '\n')
  
  # Output
  names(resultList) <- rownames(countdata)
  noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
  if (length(which(noErr)) == 0) { 
    stop("All genes returned an error. Check sufficient data in each group")
  }

  nCheat <- resultList[noErr][[1]]$predict
  outputPredict <- t(vapply(resultList[noErr], function(x) x$predict,
                            FUN.VALUE = rep(1, length(nCheat))))
  
  outLabels <- apply(modelData, 1, function(x) paste(x, collapse = "_"))
  colnames(outputPredict) <- c(paste0("y_", outLabels),
                               paste0("LCI_", outLabels),
                               paste0("UCI_", outLabels))

  if (sum(!noErr) != 0) {
    if (verbose) cat(paste0("Errors in ", sum(!noErr), " gene(s): ",
                            paste0(names(noErr)[! noErr], collapse = ", ")))
    outputErrors <- vapply(resultList[!noErr], function(x) {x$tryErrors},
                           FUN.VALUE = c("test"))
  } else {outputErrors <- c("No errors")}
  optInfo <- t(vapply(resultList[noErr], function(x) {
    setNames(x$optinfo, names(x$optinfo))
  }, FUN.VALUE = c(1, 1)))

  nCheat <- resultList[noErr][[1]]$stats
  s <- data.frame(t(vapply(resultList[noErr], function(x) {x$stats},
                           FUN.VALUE = rep(1, length(nCheat)))), 
                  check.names = FALSE)
  if(length(dispersion) > 1) s$Dispersion <- dispersion[rownames(s)]

  # Create GlmmSeq object with results
  new("GlmmSeq",
      formula = fullFormula,
      stats = s,
      predict = outputPredict,
      reducedFormula = reducedFormula,
      countdata = countdata,
      metadata = subsetMetadata,
      modelData = modelData,
      optInfo = optInfo,
      errors = outputErrors,
      variables = id
  )
}


#' Fit a glmer model for an individual gene
#'
#' @param geneList List with gene expression and optional dispersion
#' @param fullFormula the model formula. For more information of formula
#' structure see \code{\link[lme4:glmer]{lme4::glmer}}
#' @param modelData Expanded design matrix
#' @param data The sample data or metadata.
#' @param designMatrix The design matrix
#' @param control the glmer control (default = glmerControl(optimizer = 
#' "bobyqa")). For more information see
#' \code{\link[lme4:glmerControl]{lme4::glmerControl}}.
#' @param offset this can be used to specify an a priori known component to be
#'  included in the linear predictor during fitting. For more information see
#'  \code{\link[lme4:glmer]{lme4::glmer()}}.
#' @param ... Other parameters to pass to
#' \code{\link[lme4:glmer]{lme4::glmer()}}
#' @return Returns a GlmmSeq object with results for gene-wise general linear
#' mixed models
#' @importFrom MASS negative.binomial
#' @importFrom lme4 glmer fixef isSingular
#' @importFrom stats update.formula model.matrix predict setNames
#' @importFrom car Anova
#' @importFrom stats AIC complete.cases logLik reshape terms vcov predict
#' @keywords internal
#' @export
glmerApply <- function(geneList,
                       fullFormula,
                       data,
                       control,
                       modelData,
                       designMatrix,
                       offset,
                       ...) {
  data[, "count"] <- as.numeric(geneList$y)
  
  fit <- try(suppressMessages(
    lme4::glmer(fullFormula, data = data, control = control, offset = offset,
                family = geneList$family, ...)), silent = TRUE)
  
  if (class(fit) != "try-error") {
    # intercept dropped genes
    if (length(attr(fit@pp$X, "msgRankdrop")) > 0)  {
      return( list(stats = NA, predict = NA, optinfo = NA,
                   tryErrors = attr(fit@pp$X, "msgRankdrop")) )
    }
    stats <- setNames(c(AIC(fit),as.numeric(logLik(fit))),
                      c("AIC", "logLik"))
    fixedEffects <- lme4::fixef(fit)
    names(fixedEffects) <- paste0(names(fixedEffects), "_Beta")
    wald <- car::Anova(fit)
    waldtest <- setNames(c(wald[, "Chisq"], wald[, "Pr(>Chisq)"]),
                         c(paste0("Chisq_", rownames(wald)),
                           paste0("P_", rownames(wald))))
    newY <- predict(fit, newdata = modelData, re.form = NA)
    a <- designMatrix %*% vcov(fit)
    b <- as.matrix(a %*% t(designMatrix))
    predVar <- diag(b)
    newSE <- sqrt(predVar)
    newLCI <- exp(newY - newSE * 1.96)
    newUCI <- exp(newY + newSE * 1.96)
    predictdf <- c(exp(newY), newLCI, newUCI)
    singular <- as.numeric(lme4::isSingular(fit))
    conv <- length(slot(fit, "optinfo")$conv$lme4$messages)
    rm(fit, data)
    return(list(stats = c(stats, fixedEffects, waldtest),
                predict = predictdf,
                optinfo = c(singular, conv),
                tryErrors = "") )
  } else {
    return(list(stats = NA, predict = NA, optinfo = NA, tryErrors = fit[1]))
  }
}



#' Fit a glm model for an individual gene (no random effect)
#'
#' @param geneList List with gene expression and optional dispersion
#' @param fullFormula the model formula. For more information of formula
#' structure see \code{\link[stats:glm]{stats::glm}}
#' @param modelData Expanded design matrix
#' @param data The sample data or metadata.
#' @param designMatrix The design matrix
#' @param control the glm control (default = list()). For more information see
#' \code{\link[stats:glm]{stats::glm()}}.
#' @param offset this can be used to specify an a priori known component to be
#'  included in the linear predictor during fitting. For more information see
#'  \code{\link[stats:glm]{stats::glm()}}.
#' @param ... Other parameters to pass to
#' \code{\link[stats:glm]{stats::glm()}}
#' @return Returns a GlmmSeq object with results for gene-wise general linear
#' mixed models
#' @importFrom MASS negative.binomial
#' @importFrom stats glm AIC complete.cases logLik reshape terms 
#' vcov predict update.formula model.matrix predict setNames
#' @keywords internal
#' @export
glmApply <- function(geneList,
                     fullFormula,
                     data,
                     control,
                     modelData,
                     designMatrix,
                     offset,
                     ...) {
  data[, "count"] <- as.numeric(geneList$y)

  fit <- try(suppressMessages(
    glm(fullFormula, data = data, control = control, offset = offset,
        family = geneList$family, ...)), silent = TRUE)
  
  if (class(fit) != "try-error") {
    stats <- setNames(c(AIC(fit), as.numeric(logLik(fit)), fit$deviance),
                      c("AIC", "logLik", "deviance"))
    fixedEffects <- fit$coefficients
    names(fixedEffects) <- paste0(names(fixedEffects), "_Estimate")
    pValues <- summary(fit)$coefficients[, "Pr(>|t|)"]
    names(pValues) <- paste0("P_", names(pValues))
    
    newY <- predict(fit, newdata = modelData, re.form = NA)
    a <- designMatrix %*% vcov(fit)
    b <- as.matrix(a %*% t(designMatrix))
    predVar <- diag(b)
    newSE <- sqrt(predVar)
    newLCI <- exp(newY - newSE * 1.96)
    newUCI <- exp(newY + newSE * 1.96)
    predictdf <- c(exp(newY), newLCI, newUCI)
    rm(fit, data)
    return(list(stats = c(stats, fixedEffects, pValues),
                predict = predictdf,
                optinfo = c("boundary"=fit$boundary, "converged"=fit$converged),
                tryErrors = "") )
  } else {
    return(list(stats = NA, predict = NA, optinfo = NA, tryErrors = fit[1]))
  }
}
