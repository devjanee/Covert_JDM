################################################################################
# PCA Analysis for Combined Dataset (10 Contrasts)
# Creates PC1 vs PC2 and PC2 vs PC3 plots
# Style: Similar to complete_tissue_stratified_analysis_with_symbols.R
################################################################################

library(DESeq2)
library(ggplot2)
library(dplyr)

################################################################################
# STEP 1: CHECK FOR REQUIRED DATA
################################################################################

cat("=== PCA ANALYSIS FOR COMBINED DATASET ===\n\n")

if(!exists("dds_combined")) {
  stop("ERROR: dds_combined not found! Please run extract_10_contrasts.R first")
}

cat("Found DESeq2 object with", nrow(dds_combined), "genes and", 
    ncol(dds_combined), "samples\n\n")

################################################################################
# STEP 2: VARIANCE STABILIZING TRANSFORMATION
################################################################################

cat("Performing variance stabilizing transformation...\n")
vsd <- vst(dds_combined, blind = TRUE)
cat("VST complete\n\n")

################################################################################
# STEP 3: COMPUTE PCA
################################################################################

cat("Computing PCA...\n")

# Get sample information
sample_info <- as.data.frame(colData(dds_combined))

# Compute PCA manually to access all PCs
pca_obj <- prcomp(t(assay(vsd)))

# Calculate percent variance for all PCs
percentVar <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)), 1)

cat("Variance explained by top 5 PCs:\n")
for(i in 1:min(5, length(percentVar))) {
  cat(sprintf("  PC%d: %.1f%%\n", i, percentVar[i]))
}
cat("\n")

# Create PCA data frame
pca_data <- data.frame(
  Sample = rownames(pca_obj$x),
  PC1 = pca_obj$x[, 1],
  PC2 = pca_obj$x[, 2],
  PC3 = pca_obj$x[, 3],
  Group = sample_info$Group
)

cat("PCA computed successfully\n\n")

################################################################################
# STEP 4: CREATE PC1 vs PC2 PLOT
################################################################################

cat("Creating PC1 vs PC2 plot...\n")

# Define colors for groups
group_colors <- c(
  "Muscle-JDM" = "#C62828",
  "Muscle-Healthy" = "#2E7D32",
  "Myobundle-Control" = "#757575",
  "Myobundle-INFa" = "#1976D2",
  "Myobundle-INFb" = "#F57C00"
)

# Create PC1 vs PC2 plot
p_pc1_pc2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = group_colors,
                     name = "Group",
                     labels = c("Muscle-JDM" = "Muscle - JDM",
                               "Muscle-Healthy" = "Muscle - Healthy",
                               "Myobundle-Control" = "Myobundle - Control",
                               "Myobundle-INFa" = "Myobundle - IFNα",
                               "Myobundle-INFb" = "Myobundle - IFNβ")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - Combined Dataset (All Samples)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

print(p_pc1_pc2)
ggsave("PCA_Combined_PC1_vs_PC2.pdf", p_pc1_pc2, width = 10, height = 7)

cat("PC1 vs PC2 plot saved\n\n")

################################################################################
# STEP 5: CREATE PC2 vs PC3 PLOT
################################################################################

cat("Creating PC2 vs PC3 plot...\n")

# Create PC2 vs PC3 plot
p_pc2_pc3 <- ggplot(pca_data, aes(x = PC2, y = PC3, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = group_colors,
                     name = "Group",
                     labels = c("Muscle-JDM" = "Muscle - JDM",
                               "Muscle-Healthy" = "Muscle - Healthy",
                               "Myobundle-Control" = "Myobundle - Control",
                               "Myobundle-INFa" = "Myobundle - IFNα",
                               "Myobundle-INFb" = "Myobundle - IFNβ")) +
  xlab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylab(paste0("PC3: ", percentVar[3], "% variance")) +
  ggtitle("PCA - Combined Dataset (All Samples)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

print(p_pc2_pc3)
ggsave("PCA_Combined_PC2_vs_PC3.pdf", p_pc2_pc3, width = 10, height = 7)

cat("PC2 vs PC3 plot saved\n\n")

################################################################################
# STEP 6: SAVE PCA DATA
################################################################################

cat("Saving PCA data table...\n")

# Create comprehensive PCA data table
pca_data_full <- data.frame(
  Sample = rownames(pca_obj$x),
  Group = sample_info$Group,
  PC1 = pca_obj$x[, 1],
  PC2 = pca_obj$x[, 2],
  PC3 = pca_obj$x[, 3],
  PC4 = pca_obj$x[, 4],
  PC5 = pca_obj$x[, 5]
)

write.csv(pca_data_full, "PCA_Combined_Coordinates.csv", row.names = FALSE)

cat("PCA data saved\n\n")

################################################################################
# STEP 7: VARIANCE EXPLAINED SUMMARY
################################################################################

cat("Creating variance explained summary...\n")

# Create variance explained table
variance_summary <- data.frame(
  PC = paste0("PC", 1:10),
  Variance_Percent = percentVar[1:10],
  Cumulative_Percent = cumsum(percentVar[1:10])
)

write.csv(variance_summary, "PCA_Variance_Explained.csv", row.names = FALSE)

cat("\n=== VARIANCE EXPLAINED (Top 10 PCs) ===\n")
print(variance_summary)
cat("\n")

################################################################################
# STEP 8: CREATE COMBINED PLOT (OPTIONAL)
################################################################################

cat("Creating combined PC plot panel...\n")

library(gridExtra)

# Combine both plots side by side
p_combined <- grid.arrange(p_pc1_pc2, p_pc2_pc3, ncol = 2)

ggsave("PCA_Combined_Panel.pdf", p_combined, width = 18, height = 7)

cat("Combined panel plot saved\n\n")

################################################################################
# COMPLETE
################################################################################

cat("========================================\n")
cat("PCA ANALYSIS COMPLETE!\n")
cat("========================================\n\n")

cat("Generated files:\n")
cat("1. PCA_Combined_PC1_vs_PC2.pdf\n")
cat("2. PCA_Combined_PC2_vs_PC3.pdf\n")
cat("3. PCA_Combined_Panel.pdf (both plots side-by-side)\n")
cat("4. PCA_Combined_Coordinates.csv (PC scores for all samples)\n")
cat("5. PCA_Variance_Explained.csv (variance table)\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("- PC1 explains %.1f%% of variance\n", percentVar[1]))
cat(sprintf("- PC2 explains %.1f%% of variance\n", percentVar[2]))
cat(sprintf("- PC3 explains %.1f%% of variance\n", percentVar[3]))
cat(sprintf("- Top 3 PCs explain %.1f%% of total variance\n", 
            sum(percentVar[1:3])))
