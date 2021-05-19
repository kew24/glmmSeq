## ----setup, include=FALSE---------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 8, fig.height = 6)
options(width=96)
library(kableExtra)

## ---- eval=FALSE------------------------------------------------------------------------------
#  install.packages("glmmSeq")

## ---- eval=FALSE------------------------------------------------------------------------------
#  devtools::install_github("KatrionaGoldmann/glmmSeq")

## ---- eval=FALSE------------------------------------------------------------------------------
#  functions = list.files("./R", full.names = TRUE)
#  invisible(lapply(functions, source))

## ---- eval=FALSE------------------------------------------------------------------------------
#  # Install CRAN packages
#  invisible(lapply(c("MASS", "car", "ggplot2", "ggpubr", "lme4", "methods",
#                     "parallel", "plotly", "stats", "gghalves"),
#                   function(p){
#                     if(! p %in% rownames(installed.packages())) {
#                       install.packages(p)
#                     }
#                     library(p, character.only=TRUE)
#                   }))
#  
#  # Install BioConductor packages
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  invisible(lapply(c("qvalue"), function(p){
#    if(! p %in% rownames(installed.packages())) BiocManager::install(p)
#    library(p, character.only=TRUE)
#  }))
#  

## ---- message=FALSE, warning=FALSE------------------------------------------------------------
library(glmmSeq)
set.seed(1234)

## ---------------------------------------------------------------------------------------------
data(PEAC_minimal_load)

## ---------------------------------------------------------------------------------------------
metadata$EULAR_binary  = NA
metadata$EULAR_binary[metadata$EULAR_6m %in%
                        c("Good responder", "Moderate responder" )] = "responder"
metadata$EULAR_binary[metadata$EULAR_6m %in% c("Non responder")] = "non_responder"
metadata = metadata[! is.na(metadata$EULAR_binary), ]

kable(head(metadata), row.names = F) %>% kable_styling()

## ---------------------------------------------------------------------------------------------
tpm = tpm[, metadata$SAMID]
kable(head(tpm)) %>% kable_styling() %>%
  scroll_box(width = "100%")

## ---------------------------------------------------------------------------------------------
disp <- apply(tpm, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
  })

head(disp)

## ---- message=FALSE---------------------------------------------------------------------------
disp  <- setNames(edgeR::estimateDisp(tpm)$tagwise.dispersion, rownames(tpm))

head(disp)

## ---- eval=FALSE------------------------------------------------------------------------------
#  dds <- DESeqDataSetFromTximport(txi = txi, colData = metadata, design = ~ 1)
#  dds <- DESeq(dds)
#  dispersions <- setNames(dispersions(dds), rownames(txi$counts))

## ---------------------------------------------------------------------------------------------
sizeFactors <- colSums(tpm)  
sizeFactors <- sizeFactors / mean(sizeFactors)  # normalise

head(sizeFactors)

## ---- eval=FALSE------------------------------------------------------------------------------
#  sizeFactors <- calcNormFactors(counts, method="TMM")

## ---- eval=FALSE------------------------------------------------------------------------------
#  sizeFactors <- estimateSizeFactorsForMatrix(counts)

## ---- warning=FALSE---------------------------------------------------------------------------
results <- glmmSeq(modelFormula = ~ Timepoint * EULAR_6m + (1 | PATID),
                  id = "PATID",
                  countdata = tpm,
                  metadata = metadata,
                  dispersion = disp,
                  family = NULL,
                  control=list(),
                  removeDuplicatedMeasures = FALSE,
                  removeSingles=FALSE,
                  progress=TRUE,
                  cores = 1)

