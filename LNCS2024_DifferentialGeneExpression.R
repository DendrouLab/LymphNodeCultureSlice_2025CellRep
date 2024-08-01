# Updated 31-07-2024
# Differential abundance (DA) analysis
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

#### Differential abundance analysis between whole lymph node and control slice lymph node
# Import data
whole <- readRDS("whole_sce.rds")
slice <- readRDS("slice_sce.rds")

# Subsetting the control samples into a slice_ctrl object
slice_ctrl <- slice[, slice$sample_type  %in% c("sCTRL")]
slice_ctrl$sampletype <- droplevels(slice_ctrl$sample_type)
slice_ctrl$sample_id <- droplevels(slice_ctrl$sample_id)

# Excluding T, B and plasma cells prior to DA analysis
whole_da <- whole[,whole$celltypes %in% c("CD16hi_NK", "CD16low_NK", "DC","Endothelial","SELL+ID2hi_ILC","Mast_cells","Monocyte_Macrophages", "ILC3", "HLADRhi_ILC3", "Neutrophil_Granulocytes","Proliferating_cells","Stromal","pDC")]
slice_da <- slice_ctrl[,slice_ctrl$celltypes %in% c("CD25+IL22+_ILC3", "CD16hi_NK", "CD16low_NK", "DC", "Endothelial", "HSPhi_ILC", "SELL+ID2hi_ILC", "Mast_cells", "Monocyte_Macrophages", "ILC3", "HLADRhi_ILC3", "CD25+CXCL8+_ILC", "Proliferating_cells", "Stromal", "pDC")]

whole_da$celltypes <- droplevels(whole_da$celltypes)
slice_da$celltypes <- droplevels(slice_da$celltypes)

# Calculating counts per sample
whole_da <- colData(whole_da)
slice_da <- colData(slice_da)
diff_abundance_w <- as.data.frame(table(whole_da$sample_id, whole_da$celltypes))
diff_abundance_s <- as.data.frame(table(slice_da$sample_id, slice_da$celltypes))
diff_abundance <- rbind(diff_abundance_w, diff_abundance_s)
colnames(diff_abundance) <- c("sample_id", "cell_type", "number_of_cells")
diff_abundance <- separate(data = diff_abundance, col = "sample_id", c('donor_id', 'sample_type'), sep = "_", remove = FALSE)
diff_abundance$donor_id <- as.factor(diff_abundance$donor_id)
diff_abundance$sample_id <- as.factor(diff_abundance$sample_id)
diff_abundance$sample_type <- as.factor(diff_abundance$sample_type)

# Adding in a column with the total number of cells per sample
cells_per_sample <- rbind(data.frame(table(slice_da$sample_id)), data.frame(table(whole_da$sample_id)))
colnames(cells_per_sample) <- c("sample_id", "cells_per_sample")
matching_rows <- match(diff_abundance$sample_id, cells_per_sample$sample_id)
diff_abundance$cells_per_sample <- cells_per_sample$cells_per_sample[matching_rows]

# Calculating the proportion of each cell type per sample
diff_abundance$percentage_cells <- diff_abundance$number_of_cells/diff_abundance$cells_per_sample
filtered_diff_abundance <- diff_abundance %>% filter(!is.na(percentage_cells))

# Subsetting the table to only include those cell types that are common between the whole LN and slice culture datasets
diff_abundance1 <- subset(diff_abundance, cell_type %in% c("Stromal", "Endothelial", "CD16hi_NK", "CD16low_NK", "SELL+ID2hi_ILC", "ILC3", "HLADRhi_ILC3", "DC", "pDC", "Monocyte_Macrophages", "Mast_cells", "Proliferating_cells"))

# Quantifying the significance by performing an anova test where the overall error sum of squares and degrees of freedom are calculated from a linear model
cd19.model <- lm(percentage_cells ~ cell_type + sample_type, data = diff_abundance1)
significance_table <- diff_abundance1 %>%
  group_by(cell_type) %>%
  anova_test(percentage_cells ~ sample_type, error = cd19.model)

# Ploting the per sample per cell type abundances on a barplot
filtered_diff_abundance1 <- filtered_diff_abundance
filtered_diff_abundance1$sample_type <- factor(filtered_diff_abundance1$sample_type, 
                                              levels = c("wCTRL", "sCTRL"))
ggplot(filtered_diff_abundance1, aes(x = factor(cell_type, level = c("Stromal", "Endothelial", "CD16hi_NK", "CD16low_NK", "SELL+ID2hi_ILC", "ILC3", "HLADRhi_ILC3", "CD25+IL22+_ILC3", "CD25+CXCL8+_ILC", "HSPhi_ILC", "DC", "pDC", "Monocyte_Macrophages", "Neutrophil_Granulocytes", "Mast_cells", "Proliferating_cells")), 
                                    y = percentage_cells)) +
  geom_boxplot(aes(colour = sample_type, fill = sample_type), outlier.shape = NA) + 
  scale_color_manual(values = c("#000000", "#000000"))+
  scale_fill_manual(values = c("white", "#E75480")) +
  geom_point(aes(color = donor_id, group = combined_type), size = 1, position = position_dodge(width = 0.75), colour = "black") +
  labs(x = "Cell Type", y = "Proportion of cells (average across donors)", colour = "Sample Type") +
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(angle = 55, hjust = 1),
        panel.background = element_rect(fill = "white", colour = "darkgray", size = 2),
        panel.grid.major = element_line(colour = "white", size = 0.3)) 
