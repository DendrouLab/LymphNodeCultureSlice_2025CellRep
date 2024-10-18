# Updated 17-10-2024
# Differential gene expression (DA) analysis
# Written by Nitya Gupta and Jacqueline HY Siu

# Libraries
library(scater)
library(SingleCellExperiment)
library(scran)
library(scuttle)
library(speckle)
library(limma)
library(ggplot2)
library(dittoSeq)
library(ggplot2)
library(ComplexHeatmap)
library(SingleR)
library(ggalluvial)
library(circlize)
library(RColorBrewer)
library(rstatix)
library(tidyr)
library(dplyr)
library(emmeans)

## Differential Gene Expression
# import data
sce <- readRDS("file.rds")

# conducting a pairwise t test on the normalised log genes expression values between the CTRL and LMQ groupd for each cell type 
sce$clusterid.sampletype <- paste(sce$celltypes, sce$sample_type, sep = "_")
out <- pairwiseTTests(logcounts(sce), groups=sce$clusterid.sampletype, block = sce$celltypes)
pairings <- out$pairs
pairings$combined <- paste0(pairings$first, "-", pairings$second)
pairings$rowid <- 1:nrow(pairings)
temp <- {}

for(i in as.factor(unique(sce$celltypes))){
  temp <- c(temp, paste0(i, "_sLMQ-", i, "_sCTRL"))
}

keep <- subset(pairings, combined %in% temp)
stats_keep <- out[["statistics"]][keep$rowid]
pairs_keep <- out[["pairs"]][keep$rowid,]
cbm_all <- scran::combineMarkers(
  de.lists = stats_keep, pairs = pairs_keep)

#Save CSV and VolcanoPlot
temp_names <- names(cbm_all)
temp_names <- gsub("_sLMQ", "", temp_names)
for (i in seq(cbm_all)){
  temp <- cbm_all[[i]]
  write.csv(temp, file = file.path("SCE_DEG", paste0("DEGscran_cluster",temp_names[i], ".csv")))
  EnhancedVolcano(temp, lab = rownames(temp), x = 'summary.logFC', y = 'p.value', FCcutoff = 0.25, drawConnectors = TRUE, title = paste0(temp_names[i]), subtitle = 'sCTRL (LH), sLMQ (RH)')
  ggsave(file.path(out_dir, "SCE_DEG_res06", paste0("SCEDEG_VolcanoPlot", temp_names[i],".png")), width = 8, height = 8)
}
