#' Glmm for sequencing results of a single gene
#'
#' @param modelFormula the model formula. For more information of formula
#' structure see \code{\link[lme4:glmer]{lme4::glmer}}.
#' @param countdata the sequencing data
#' @param gene the row name in countdata to be used
#' @param metadata a data frame of sample information
#' @param id Column name in metadata which contains the sample IDs to be used
#' in pairing
#' @param dispersion a numeric for the gene dispersion
#' @param sizeFactors size factors (default=NULL). If provided the glmer offset
#' is set to log(sizeFactors).For more information see
#'  \code{\link[lme4:glmer]{lme4::glmer()}}
#' @param reducedFormula Reduced design formula (default="")
#' @param modelData something something
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
#' @param zeroCount numerical value to offset zeroes for the purpose of log
#' (default=0.125)
#' @param removeDuplicatedMeasures whether to remove duplicated
#' conditions/repeated measurements for a given time point (default=FALSE).
#' Only used if glmOnly = FALSE. 
#' @param removeSingles whether to remove individuals with only one measurement
#' (default=FALSE). Only used if glmOnly = FALSE. 
#' @param verbose Logical whether to display messaging (default=FALSE)
#' @param ... Other parameters to pass to
#' \code{\link[lme4:glmer]{lme4::glmer}}.
#' @return Returns the fit for the general linear mixed model of a single gene
#' @importFrom MASS negative.binomial
#' @importFrom lme4 subbars findbars glmer fixef glmerControl nobars
#' @importFrom stats update.formula model.matrix predict setNames
#' @importFrom car Anova
#' @importFrom methods slot new
#' @importFrom stats AIC complete.cases logLik reshape terms vcov
#' @export
#' @examples
#' data(PEAC_minimal_load)
#' disp <- apply(tpm, 1, function(x) {
#' (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
#' })
#' MS4A1fit <- glmmGene(~ Timepoint * EULAR_6m + (1 | PATID),
#'                       gene = "MS4A1",
#'                       id = "PATID",
#'                       countdata = tpm,
#'                       metadata = metadata,
#'                       dispersion = disp["MS4A1"],
#'                       verbose=FALSE)
#'
#' MS4A1fit

glmmGene <- function(modelFormula,
                     countdata,
                     gene,
                     metadata,
                     id,
                     dispersion,
                     sizeFactors = NULL,
                     reducedFormula = "",
                     modelData = NULL,
                     glmOnly = FALSE,
                     control = glmerControl(optimizer = "bobyqa"),
                     family = NULL,
                     zeroCount = 0.125,
                     removeDuplicatedMeasures = FALSE,
                     removeSingles = FALSE,
                     verbose = FALSE,
                     ...){
  
  # Catch errors
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
  if (! gene %in% rownames(countdata)) {
    stop("gene must be in rownames(countdata)")
  }
  if(is.null(family)){
    family <- MASS::negative.binomial(theta = 1/dispersion)
  }
  if (! is.numeric(zeroCount)) stop("zeroCount must be numeric")
  if (zeroCount < 0) stop("zeroCount must be >= 0")
  if (zeroCount > 0) countdata[countdata == 0] <- zeroCount
  
  # Manipulate formulae
  fullFormula <- update.formula(modelFormula, count ~ ., simplify = FALSE)
  nonRandomFormula <- subbars(modelFormula)
  variables <- rownames(attr(terms(nonRandomFormula), "factors"))
  subsetMetadata <- metadata[, variables]
  
  if(! glmOnly) {
    if(! id %in% colnames(metadata)) 
      stop("id must be a column name in metadata")
    ids <- as.character(metadata[, id])
    
    
    # Option to subset to remove duplicated timepoints
    if (removeDuplicatedMeasures) {
      # Check the distribution for duplicates
      check <- data.frame(table(droplevels(subsetMetadata)))
      check <- check[! check$Freq %in% c(0, 1), ]
      if (nrow(check) > 0) {
        mCheck <- as.character(
          apply(subsetMetadata[, variables], 1, function(x) {
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
      stop("Alignment error")
    }
    if (verbose) cat(paste0("\nn=", length(ids), " samples, ",
                            length(unique(ids)), " individuals\n"))
  }
  if (!is.null(sizeFactors)) offset <- log(sizeFactors) else offset <- NULL
  
  # setup model prediction
  if (reducedFormula == "") reducedFormula <- nobars(modelFormula)
  if (is.null(modelData)) {
    reducedVars <- rownames(attr(terms(reducedFormula), "factors"))
    varLevels <- lapply(reducedVars, function(x) {
      if (class(metadata[, x]) == "factor") return(levels(metadata[, x]))
      sort(unique(metadata[, x]))
    })
    modelData <- expand.grid(varLevels)
    colnames(modelData) <- reducedVars
  }
  
  data <- subsetMetadata
  data[, "count"] <- as.numeric(countdata[gene, ])
  
  
  # Perform a fit
  if(glmOnly) {
    functionUsed <- glmApply 
    if("glmerControl" %in% class(control)) {
      control <- list()
      message("glmer control passed into glm, using control <- list() instead")
    }
    
    fit <- try(glm(fullFormula, data = data, control = control, 
                   offset = offset, family = family, ...), silent = TRUE)
    
  } else {
    
    if(! "glmerControl" %in% class(control)) {
      control <- glmerControl()
      message(paste("Control not of class 'glmer control', using control <-", 
                    "glmerControl() instead"))
    }
    fit <- try(glmer(fullFormula, data = data, control = control, 
                     offset = offset, family=family, ...), silent=FALSE)
  }
  
  
  return(fit)
  
}
