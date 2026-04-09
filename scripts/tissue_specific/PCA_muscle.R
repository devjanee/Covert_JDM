cat("Computing PCA...\n")

# Get sample information
sample_info <- as.data.frame(colData(dds_muscle))
sample_info$Site <- ifelse(sample_info$Location == "Duke", "Site 1", "Site 2")
# Compute PCA manually to access all PCs
pca_obj <- prcomp(t(assay(vsd_muscle)))

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
  Condition = sample_info$Condition,
  Site = sample_info$Site
)

cat("PCA computed successfully\n\n")

################################################################################
# STEP 4: CREATE PC1 vs PC2 PLOT
################################################################################

cat("Creating PC1 vs PC2 plot...\n")

# Define colors for groups
group_colors <- c(
  "JDMS" = "#FF61CC",
  "Healthy" = "#00B8E7"
)

# Create PC1 vs PC2 plot
p_pc1_pc2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Site)) +
  geom_point(size = 6, alpha = 0.8) +
  scale_color_manual(values = group_colors,
                     name = "Condition",
                     labels = c("JDMS" = "JDM",
                                "Healthy" = "Healthy"
                                )) +
  scale_shape_manual(values = c("Site 1" = 16, "Site 2" = 17),
                     name = "Location",
                     labels = c("Site 1" = "Site 1",
                                "Site 2" = "Site 2")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 30) +
  theme(legend.position = "right")

print(p_pc1_pc2)
ggsave("~/Documents/Covert_11358/PCA/muscle/20250202_PCA_muscle_PC1_vs_PC2.pdf", p_pc1_pc2, width = 10, height = 7)

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
