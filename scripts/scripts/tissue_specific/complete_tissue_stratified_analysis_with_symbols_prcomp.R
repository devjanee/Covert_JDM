################################################################################
# Complete Tissue-Stratified Analysis with Gene Symbols
# Part 1: Patient Muscle Analysis (Duke + UCSF combined)
# Part 2: Myobundle Analysis (Duke only, separate normalization)
# Part 3: Overlap comparison
# UPDATED: Gene symbols added to all results and gene lists
################################################################################

library(DESeq2)
library(sva)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(dplyr)
library(RColorBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)

################################################################################
# CRITICAL POINT:
# DESeq2 normalization is calculated across ALL samples in the input matrix
# Therefore, patient muscle and myobundles MUST be analyzed separately
# to avoid cross-contamination of normalization
################################################################################

cat("=== TISSUE-STRATIFIED ANALYSIS (SEPARATE NORMALIZATIONS) ===\n\n")

################################################################################
# FUNCTION: ADD GENE SYMBOLS TO RESULTS
################################################################################

add_gene_symbols <- function(results_df, gene_col = "gene_id") {
  cat("Adding gene symbols...\n")
  
  # Get ENSEMBL IDs
  ensembl_ids <- results_df[[gene_col]]
  
  # Convert to gene symbols
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  
  # Add gene symbols - keep ENSEMBL if no symbol available
  results_df$gene_symbol <- ifelse(is.na(gene_symbols), 
                                   results_df[[gene_col]], 
                                   gene_symbols)
  
  # Reorder columns to put gene_symbol right after gene_id
  col_order <- c(gene_col, "gene_symbol", 
                 setdiff(names(results_df), c(gene_col, "gene_symbol")))
  results_df <- results_df[, col_order]
  
  return(results_df)
}

################################################################################
# 1. LOAD DATA
################################################################################

cat("Loading data...\n")

sample_sheet <- read.csv("~/Documents/Covert_11358/Covert_Sample ID list_10.10.25_DSL.csv",
                         check.names = FALSE)
counts <- read.csv("~/Documents/Covert_11358/20250605_counts_filtered.csv",
                   row.names = 1, check.names = FALSE)

# Match samples
common_samples <- intersect(colnames(counts), sample_sheet$`Sample name for RNAseq`)
counts <- counts[, common_samples]
sample_sheet <- sample_sheet[sample_sheet$`Sample name for RNAseq` %in% common_samples, ]
sample_sheet <- sample_sheet[match(colnames(counts), sample_sheet$`Sample name for RNAseq`), ]

cat("Total samples loaded:", ncol(counts), "\n")
cat("Total genes:", nrow(counts), "\n\n")

################################################################################
# PART 1: PATIENT MUSCLE ANALYSIS (DUKE + UCSF)
################################################################################

cat("============================================================\n")
cat("PART 1: PATIENT MUSCLE ANALYSIS (DUKE + UCSF COMBINED)\n")
cat("============================================================\n\n")

# Subset to muscle samples only
muscle_idx <- sample_sheet$Tissue == "Muscle"
muscle_counts <- counts[, muscle_idx]
muscle_sheet <- sample_sheet[muscle_idx, ]

cat("Patient muscle samples:\n")
cat("  Total:", ncol(muscle_counts), "\n")
cat("  Duke:", sum(muscle_sheet$Location == "Duke"), "\n")
cat("  UCSF:", sum(muscle_sheet$Location == "UCSF"), "\n")
cat("  Healthy:", sum(muscle_sheet$Condition == "Healthy"), "\n")
cat("  JDMS:", sum(muscle_sheet$Condition == "JDMS"), "\n\n")

cat("Sample distribution:\n")
print(table(muscle_sheet$Condition, muscle_sheet$Location))
cat("\n")

# Apply ComBat-seq to remove site batch effect
cat("Applying ComBat-seq to remove site batch effect...\n")

muscle_counts_corrected <- ComBat_seq(
  counts = as.matrix(muscle_counts),
  batch = muscle_sheet$Location,
  group = muscle_sheet$Condition
)

cat("ComBat-seq complete\n\n")

# Create DESeq2 object for patient muscle
cat("Creating DESeq2 object for patient muscle...\n")

muscle_sheet$Condition <- factor(muscle_sheet$Condition, levels = c("Healthy", "JDMS"))

dds_muscle <- DESeqDataSetFromMatrix(
  countData = muscle_counts_corrected,
  colData = muscle_sheet,
  design = ~ Condition
)

cat("DESeq2 object dimensions:", nrow(dds_muscle), "genes x", ncol(dds_muscle), "samples\n\n")

# Pre-filter low count genes
cat("Filtering low-count genes...\n")
keep_muscle <- rowSums(counts(dds_muscle)) >= 10
dds_muscle <- dds_muscle[keep_muscle, ]
cat("Retained", nrow(dds_muscle), "genes after filtering\n\n")

# Run DESeq2
cat("Running DESeq2 on patient muscle...\n")
dds_muscle <- DESeq(dds_muscle)
cat("DESeq2 complete!\n\n")

# Check results names
cat("Available contrasts:\n")
print(resultsNames(dds_muscle))
cat("\n")

# Extract results: JDMS vs Healthy
res_muscle_jdms <- results(dds_muscle,
                           contrast = c("Condition", "JDMS", "Healthy"),
                           alpha = 0.05)

cat("PATIENT MUSCLE: JDMS vs Healthy\n")
print(summary(res_muscle_jdms))
cat("\n")

# Count significant genes
sig_muscle_total <- sum(res_muscle_jdms$padj < 0.05, na.rm = TRUE)
sig_muscle_up <- sum(res_muscle_jdms$padj < 0.05 & res_muscle_jdms$log2FoldChange > 0, na.rm = TRUE)
sig_muscle_down <- sum(res_muscle_jdms$padj < 0.05 & res_muscle_jdms$log2FoldChange < 0, na.rm = TRUE)

cat("Summary of significant genes (padj < 0.05):\n")
cat("  Total:", sig_muscle_total, "\n")
cat("  Upregulated:", sig_muscle_up, "\n")
cat("  Downregulated:", sig_muscle_down, "\n\n")

# Convert to data frame
res_muscle_jdms_df <- as.data.frame(res_muscle_jdms)
res_muscle_jdms_df$gene_id <- rownames(res_muscle_jdms_df)
res_muscle_jdms_df <- res_muscle_jdms_df[, c("gene_id", "baseMean", "log2FoldChange", 
                                              "lfcSE", "stat", "pvalue", "padj")]

# ADD GENE SYMBOLS
res_muscle_jdms_df <- add_gene_symbols(res_muscle_jdms_df)

# Save results
write.csv(res_muscle_jdms_df, 
          "DESeq2_PatientMuscle_JDMS_vs_Healthy.csv",
          row.names = FALSE)

cat("Patient muscle results saved (with gene symbols)\n\n")

# PCA for patient muscle
cat("Generating PCA for patient muscle...\n")
vsd_muscle <- vst(dds_muscle, blind = TRUE)

pca_data <- data.frame(
  Sample = rownames(pca_obj$x),
  PC1 = pca_obj$x[, 1],
  PC2 = pca_obj$x[, 2],
  PC3 = pca_obj$x[, 3],
  Group = sample_info$Group
)
pcaData_muscle <- plotPCA(vsd_muscle, 
                          intgroup = c("Condition", "Location"), 
                          returnData = TRUE)
percentVar_muscle <- round(100 * attr(pcaData_muscle, "percentVar"))

# Calculate batch effect in PC1
pca_obj_muscle <- prcomp(t(assay(vsd_muscle)))
pc1_scores_muscle <- pca_obj_muscle$x[, 1]
batch_r2_muscle <- summary(lm(pc1_scores_muscle ~ muscle_sheet$Location))$r.squared * 100

cat("PC1 variance explained by Location:", round(batch_r2_muscle, 1), "%\n\n")

p_muscle <- ggplot(pcaData_muscle, aes(PC1, PC2, color = Condition, shape = Location)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Healthy" = "#2E7D32", "JDMS" = "#C62828")) +
  scale_shape_manual(values = c("Duke" = 16, "UCSF" = 17)) +
  xlab(paste0("PC1: ", percentVar_muscle[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_muscle[2], "% variance")) +
  ggtitle("PCA - Patient Muscle (Duke + UCSF, ComBat-corrected)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_muscle)
ggsave("PCA_PatientMuscle_Combined.pdf", p_muscle, width = 8, height = 6)

cat("Patient muscle PCA saved\n\n")

# Sample distance heatmap
cat("Creating sample distance heatmap for patient muscle...\n")

sampleDists_muscle <- dist(t(assay(vsd_muscle)))
sampleDistMatrix_muscle <- as.matrix(sampleDists_muscle)

annotation_col_muscle <- data.frame(
  Condition = muscle_sheet$Condition,
  Location = muscle_sheet$Location,
  row.names = colnames(vsd_muscle)
)

ann_colors_muscle <- list(
  Condition = c("Healthy" = "#2E7D32", "JDMS" = "#C62828"),
  Location = c("Duke" = "#1976D2", "UCSF" = "#F57C00")
)

pdf("Heatmap_PatientMuscle_SampleDistances.pdf", width = 10, height = 9)
pheatmap(sampleDistMatrix_muscle,
         clustering_distance_rows = sampleDists_muscle,
         clustering_distance_cols = sampleDists_muscle,
         annotation_col = annotation_col_muscle,
         annotation_colors = ann_colors_muscle,
         main = "Sample Distances - Patient Muscle")
dev.off()

cat("Sample distance heatmap saved\n\n")

# Volcano plot
cat("Creating volcano plot for patient muscle...\n")

res_muscle_plot <- res_muscle_jdms_df
res_muscle_plot$significant <- ifelse(res_muscle_plot$padj < 0.05 & 
                                      abs(res_muscle_plot$log2FoldChange) > 1, 
                                      "Significant", "Not Significant")
res_muscle_plot$significant[is.na(res_muscle_plot$significant)] <- "Not Significant"

p_volcano_muscle <- ggplot(res_muscle_plot, 
                           aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(title = "Patient Muscle: JDMS vs Healthy",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Volcano_PatientMuscle_JDMS_vs_Healthy.pdf", p_volcano_muscle, width = 8, height = 6)

cat("Volcano plot saved\n\n")

################################################################################
# PART 2: MYOBUNDLE ANALYSIS (DUKE ONLY, SEPARATE NORMALIZATION)
################################################################################

cat("============================================================\n")
cat("PART 2: MYOBUNDLE ANALYSIS (DUKE ONLY)\n")
cat("============================================================\n\n")

# Subset to myobundles only
myobundle_idx <- sample_sheet$Tissue == "Myobundle"
myobundle_counts <- counts[, myobundle_idx]
myobundle_sheet <- sample_sheet[myobundle_idx, ]

cat("Myobundle samples:\n")
cat("  Total:", ncol(myobundle_counts), "\n")
cat("  Location:", unique(myobundle_sheet$Location), "\n")
cat("  NT:", sum(myobundle_sheet$Condition == "NT"), "\n")
cat("  IFNa (all doses):", sum(grepl("IFNa", myobundle_sheet$Condition)), "\n")
cat("  IFNb (all doses):", sum(grepl("IFNb", myobundle_sheet$Condition)), "\n\n")

cat("Sample distribution by condition:\n")
print(table(myobundle_sheet$Condition))
cat("\n")

cat("Sample distribution by donor:\n")
print(table(myobundle_sheet$Individual, myobundle_sheet$Condition))
cat("\n")

# Create grouping variable (lump IFN doses together)
myobundle_sheet$group <- NA
myobundle_sheet$group[myobundle_sheet$Condition == "NT"] <- "NT"
myobundle_sheet$group[myobundle_sheet$Condition %in% c("IFNa5", "IFNa10", "IFNa20")] <- "IFNa"
myobundle_sheet$group[myobundle_sheet$Condition %in% c("IFNb5", "IFNb10", "IFNb20")] <- "IFNb"

myobundle_sheet$group <- factor(myobundle_sheet$group, levels = c("NT", "IFNa", "IFNb"))

# Also keep detailed groups for PCA
myobundle_sheet$group_detailed <- factor(myobundle_sheet$Condition)

cat("Myobundle groups (lumped):\n")
print(table(myobundle_sheet$group))
cat("\n")

# Create DESeq2 object for myobundles
cat("Creating DESeq2 object for myobundles...\n")

dds_myobundle <- DESeqDataSetFromMatrix(
  countData = myobundle_counts,
  colData = myobundle_sheet,
  design = ~ group
)

cat("DESeq2 object dimensions:", nrow(dds_myobundle), "genes x", ncol(dds_myobundle), "samples\n\n")

# Pre-filter low count genes
cat("Filtering low-count genes...\n")
keep_myobundle <- rowSums(counts(dds_myobundle)) >= 10
dds_myobundle <- dds_myobundle[keep_myobundle, ]
cat("Retained", nrow(dds_myobundle), "genes after filtering\n\n")

# Run DESeq2
cat("Running DESeq2 on myobundles...\n")
dds_myobundle <- DESeq(dds_myobundle)
cat("DESeq2 complete!\n\n")

# Check results names
cat("Available contrasts:\n")
print(resultsNames(dds_myobundle))
cat("\n")

# Extract results: IFNa vs NT
res_ifna_vs_nt <- results(dds_myobundle,
                          contrast = c("group", "IFNa", "NT"),
                          alpha = 0.05)

cat("MYOBUNDLES: IFNa vs NT\n")
print(summary(res_ifna_vs_nt))
cat("\n")

# Extract results: IFNb vs NT
res_ifnb_vs_nt <- results(dds_myobundle,
                          contrast = c("group", "IFNb", "NT"),
                          alpha = 0.05)

cat("MYOBUNDLES: IFNb vs NT\n")
print(summary(res_ifnb_vs_nt))
cat("\n")

# Extract results: IFNb vs IFNa (bonus comparison)
res_ifnb_vs_ifna <- results(dds_myobundle,
                            contrast = c("group", "IFNb", "IFNa"),
                            alpha = 0.05)

cat("MYOBUNDLES: IFNb vs IFNa\n")
print(summary(res_ifnb_vs_ifna))
cat("\n")

# Count significant genes
sig_ifna_total <- sum(res_ifna_vs_nt$padj < 0.05, na.rm = TRUE)
sig_ifna_up <- sum(res_ifna_vs_nt$padj < 0.05 & res_ifna_vs_nt$log2FoldChange > 0, na.rm = TRUE)
sig_ifna_down <- sum(res_ifna_vs_nt$padj < 0.05 & res_ifna_vs_nt$log2FoldChange < 0, na.rm = TRUE)

sig_ifnb_total <- sum(res_ifnb_vs_nt$padj < 0.05, na.rm = TRUE)
sig_ifnb_up <- sum(res_ifnb_vs_nt$padj < 0.05 & res_ifnb_vs_nt$log2FoldChange > 0, na.rm = TRUE)
sig_ifnb_down <- sum(res_ifnb_vs_nt$padj < 0.05 & res_ifnb_vs_nt$log2FoldChange < 0, na.rm = TRUE)

cat("Summary of significant genes:\n")
cat("IFNa vs NT:\n")
cat("  Total:", sig_ifna_total, "\n")
cat("  Upregulated:", sig_ifna_up, "\n")
cat("  Downregulated:", sig_ifna_down, "\n\n")

cat("IFNb vs NT:\n")
cat("  Total:", sig_ifnb_total, "\n")
cat("  Upregulated:", sig_ifnb_up, "\n")
cat("  Downregulated:", sig_ifnb_down, "\n\n")

# Convert to data frames
res_ifna_df <- as.data.frame(res_ifna_vs_nt)
res_ifna_df$gene_id <- rownames(res_ifna_df)
res_ifna_df <- res_ifna_df[, c("gene_id", "baseMean", "log2FoldChange", 
                               "lfcSE", "stat", "pvalue", "padj")]

res_ifnb_df <- as.data.frame(res_ifnb_vs_nt)
res_ifnb_df$gene_id <- rownames(res_ifnb_df)
res_ifnb_df <- res_ifnb_df[, c("gene_id", "baseMean", "log2FoldChange", 
                               "lfcSE", "stat", "pvalue", "padj")]

res_ifnb_vs_ifna_df <- as.data.frame(res_ifnb_vs_ifna)
res_ifnb_vs_ifna_df$gene_id <- rownames(res_ifnb_vs_ifna_df)
res_ifnb_vs_ifna_df <- res_ifnb_vs_ifna_df[, c("gene_id", "baseMean", "log2FoldChange", 
                                                "lfcSE", "stat", "pvalue", "padj")]

# ADD GENE SYMBOLS
res_ifna_df <- add_gene_symbols(res_ifna_df)
res_ifnb_df <- add_gene_symbols(res_ifnb_df)
res_ifnb_vs_ifna_df <- add_gene_symbols(res_ifnb_vs_ifna_df)

# Save results
write.csv(res_ifna_df, 
          "DESeq2_Myobundles_IFNa_vs_NT.csv",
          row.names = FALSE)
write.csv(res_ifnb_df, 
          "DESeq2_Myobundles_IFNb_vs_NT.csv",
          row.names = FALSE)
write.csv(res_ifnb_vs_ifna_df,
          "DESeq2_Myobundles_IFNb_vs_IFNa.csv",
          row.names = FALSE)

cat("Myobundle results saved (with gene symbols)\n\n")

# PCA for myobundles (with dose information)
cat("Generating PCA for myobundles...\n")
vsd_myobundle <- vst(dds_myobundle, blind = TRUE)

pcaData_myobundle <- plotPCA(vsd_myobundle, 
                             intgroup = c("group_detailed"), 
                             returnData = TRUE)
percentVar_myobundle <- round(100 * attr(pcaData_myobundle, "percentVar"))

# Color palette for all conditions
color_palette_myobundle <- c(
  "NT" = "#757575",
  "IFNa5" = "#64B5F6",
  "IFNa10" = "#1976D2",
  "IFNa20" = "#0D47A1",
  "IFNb5" = "#FFB74D",
  "IFNb10" = "#F57C00",
  "IFNb20" = "#E65100"
)

p_myobundle <- ggplot(pcaData_myobundle, aes(PC1, PC2, color = group_detailed)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = color_palette_myobundle, name = "Treatment") +
  xlab(paste0("PC1: ", percentVar_myobundle[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_myobundle[2], "% variance")) +
  ggtitle("PCA - Myobundles (Duke)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_myobundle)
ggsave("PCA_Myobundles_AllDoses.pdf", p_myobundle, width = 9, height = 6)

cat("Myobundle PCA saved\n\n")

# Sample distance heatmap
cat("Creating sample distance heatmap for myobundles...\n")

sampleDists_myobundle <- dist(t(assay(vsd_myobundle)))
sampleDistMatrix_myobundle <- as.matrix(sampleDists_myobundle)

annotation_col_myobundle <- data.frame(
  Group = myobundle_sheet$group,
  Treatment = myobundle_sheet$Condition,
  Donor = myobundle_sheet$Individual,
  row.names = colnames(vsd_myobundle)
)

ann_colors_myobundle <- list(
  Group = c("NT" = "#757575", "IFNa" = "#1976D2", "IFNb" = "#F57C00"),
  Treatment = color_palette_myobundle
)

pdf("Heatmap_Myobundles_SampleDistances.pdf", width = 11, height = 10)
pheatmap(sampleDistMatrix_myobundle,
         clustering_distance_rows = sampleDists_myobundle,
         clustering_distance_cols = sampleDists_myobundle,
         annotation_col = annotation_col_myobundle,
         annotation_colors = ann_colors_myobundle,
         main = "Sample Distances - Myobundles")
dev.off()

cat("Sample distance heatmap saved\n\n")

# Volcano plots
cat("Creating volcano plots for myobundles...\n")

# IFNa vs NT
res_ifna_plot <- res_ifna_df
res_ifna_plot$significant <- ifelse(res_ifna_plot$padj < 0.05 & 
                                    abs(res_ifna_plot$log2FoldChange) > 1, 
                                    "Significant", "Not Significant")
res_ifna_plot$significant[is.na(res_ifna_plot$significant)] <- "Not Significant"

p_volcano_ifna <- ggplot(res_ifna_plot, 
                         aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(title = "Myobundles: IFNa vs NT",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Volcano_Myobundles_IFNa_vs_NT.pdf", p_volcano_ifna, width = 8, height = 6)

# IFNb vs NT
res_ifnb_plot <- res_ifnb_df
res_ifnb_plot$significant <- ifelse(res_ifnb_plot$padj < 0.05 & 
                                    abs(res_ifnb_plot$log2FoldChange) > 1, 
                                    "Significant", "Not Significant")
res_ifnb_plot$significant[is.na(res_ifnb_plot$significant)] <- "Not Significant"

p_volcano_ifnb <- ggplot(res_ifnb_plot, 
                         aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(title = "Myobundles: IFNb vs NT",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Volcano_Myobundles_IFNb_vs_NT.pdf", p_volcano_ifnb, width = 8, height = 6)

cat("Volcano plots saved\n\n")

# Create summary statistics table
summary_stats <- data.frame(
  Comparison = c("Patient Muscle: JDMS vs Healthy",
                 "Myobundles: IFNa vs NT",
                 "Myobundles: IFNb vs NT"),
  Total_Significant = c(sig_muscle_total, sig_ifna_total, sig_ifnb_total),
  Upregulated = c(sig_muscle_up, sig_ifna_up, sig_ifnb_up),
  Downregulated = c(sig_muscle_down, sig_ifna_down, sig_ifnb_down)
)

print(summary_stats)
write.csv(summary_stats, "DESeq2_Summary_Statistics.csv", row.names = FALSE)

cat("\n")

################################################################################
# PART 3: OVERLAP ANALYSIS
################################################################################

cat("============================================================\n")
cat("PART 3: OVERLAP ANALYSIS (JDMS vs IFN)\n")
cat("============================================================\n\n")

cat("Extracting significant gene lists (padj < 0.05)...\n\n")

# JDMS signature (patient muscle) - BY GENE SYMBOL
jdms_all <- res_muscle_jdms_df$gene_symbol[
  res_muscle_jdms_df$padj < 0.05 & !is.na(res_muscle_jdms_df$padj)
]
jdms_up <- res_muscle_jdms_df$gene_symbol[
  res_muscle_jdms_df$padj < 0.05 & 
  res_muscle_jdms_df$log2FoldChange > 0 &
  !is.na(res_muscle_jdms_df$padj)
]
jdms_down <- res_muscle_jdms_df$gene_symbol[
  res_muscle_jdms_df$padj < 0.05 & 
  res_muscle_jdms_df$log2FoldChange < 0 &
  !is.na(res_muscle_jdms_df$padj)
]

# IFNa signature (myobundles) - BY GENE SYMBOL
ifna_all <- res_ifna_df$gene_symbol[
  res_ifna_df$padj < 0.05 & !is.na(res_ifna_df$padj)
]
ifna_up <- res_ifna_df$gene_symbol[
  res_ifna_df$padj < 0.05 & 
  res_ifna_df$log2FoldChange > 0 &
  !is.na(res_ifna_df$padj)
]
ifna_down <- res_ifna_df$gene_symbol[
  res_ifna_df$padj < 0.05 & 
  res_ifna_df$log2FoldChange < 0 &
  !is.na(res_ifna_df$padj)
]

# IFNb signature (myobundles) - BY GENE SYMBOL
ifnb_all <- res_ifnb_df$gene_symbol[
  res_ifnb_df$padj < 0.05 & !is.na(res_ifnb_df$padj)
]
ifnb_up <- res_ifnb_df$gene_symbol[
  res_ifnb_df$padj < 0.05 & 
  res_ifnb_df$log2FoldChange > 0 &
  !is.na(res_ifnb_df$padj)
]
ifnb_down <- res_ifnb_df$gene_symbol[
  res_ifnb_df$padj < 0.05 & 
  res_ifnb_df$log2FoldChange < 0 &
  !is.na(res_ifnb_df$padj)
]

# Print sizes
cat("GENE LIST SIZES (BY GENE SYMBOL):\n\n")
cat("JDMS SIGNATURE (Patient Muscle, Duke+UCSF):\n")
cat("  Total:", length(jdms_all), "genes\n")
cat("  Upregulated:", length(jdms_up), "genes\n")
cat("  Downregulated:", length(jdms_down), "genes\n\n")

cat("IFNa RESPONSE (Myobundles, Duke):\n")
cat("  Total:", length(ifna_all), "genes\n")
cat("  Upregulated:", length(ifna_up), "genes\n")
cat("  Downregulated:", length(ifna_down), "genes\n\n")

cat("IFNb RESPONSE (Myobundles, Duke):\n")
cat("  Total:", length(ifnb_all), "genes\n")
cat("  Upregulated:", length(ifnb_up), "genes\n")
cat("  Downregulated:", length(ifnb_down), "genes\n\n")

# Calculate overlaps
overlap_up_ifna <- intersect(jdms_up, ifna_up)
overlap_up_ifnb <- intersect(jdms_up, ifnb_up)
overlap_down_ifna <- intersect(jdms_down, ifna_down)
overlap_down_ifnb <- intersect(jdms_down, ifnb_down)
overlap_all_ifna <- intersect(jdms_all, ifna_all)
overlap_all_ifnb <- intersect(jdms_all, ifnb_all)

# Recapitulation percentages
recap_up_ifna <- length(overlap_up_ifna) / length(jdms_up) * 100
recap_up_ifnb <- length(overlap_up_ifnb) / length(jdms_up) * 100
recap_down_ifna <- length(overlap_down_ifna) / length(jdms_down) * 100
recap_down_ifnb <- length(overlap_down_ifnb) / length(jdms_down) * 100
recap_all_ifna <- length(overlap_all_ifna) / length(jdms_all) * 100
recap_all_ifnb <- length(overlap_all_ifnb) / length(jdms_all) * 100

cat("=== RECAPITULATION ANALYSIS ===\n\n")

cat("IFNa treatment recapitulates JDMS:\n")
cat(sprintf("  Upregulated: %d/%d genes (%.1f%%)\n", 
            length(overlap_up_ifna), length(jdms_up), recap_up_ifna))
cat(sprintf("  Downregulated: %d/%d genes (%.1f%%)\n", 
            length(overlap_down_ifna), length(jdms_down), recap_down_ifna))
cat(sprintf("  Total: %d/%d genes (%.1f%%)\n\n", 
            length(overlap_all_ifna), length(jdms_all), recap_all_ifna))

cat("IFNb treatment recapitulates JDMS:\n")
cat(sprintf("  Upregulated: %d/%d genes (%.1f%%)\n", 
            length(overlap_up_ifnb), length(jdms_up), recap_up_ifnb))
cat(sprintf("  Downregulated: %d/%d genes (%.1f%%)\n", 
            length(overlap_down_ifnb), length(jdms_down), recap_down_ifnb))
cat(sprintf("  Total: %d/%d genes (%.1f%%)\n\n", 
            length(overlap_all_ifnb), length(jdms_all), recap_all_ifnb))

# Which IFN is better?
if(recap_all_ifna > recap_all_ifnb) {
  cat(sprintf("CONCLUSION: IFNa better recapitulates JDMS (%.1f%% vs %.1f%%)\n\n",
              recap_all_ifna, recap_all_ifnb))
} else if(recap_all_ifnb > recap_all_ifna) {
  cat(sprintf("CONCLUSION: IFNb better recapitulates JDMS (%.1f%% vs %.1f%%)\n\n",
              recap_all_ifnb, recap_all_ifna))
} else {
  cat("CONCLUSION: IFNa and IFNb equally recapitulate JDMS\n\n")
}

# Create overlap summary table
overlap_summary <- data.frame(
  Metric = c("JDMS genes (total)",
             "JDMS genes (up)",
             "JDMS genes (down)",
             "IFNa genes (total)",
             "IFNa genes (up)",
             "IFNa genes (down)",
             "IFNb genes (total)",
             "IFNb genes (up)",
             "IFNb genes (down)",
             "JDMS-IFNa overlap (total)",
             "JDMS-IFNa overlap (up)",
             "JDMS-IFNa overlap (down)",
             "JDMS-IFNb overlap (total)",
             "JDMS-IFNb overlap (up)",
             "JDMS-IFNb overlap (down)",
             "IFNa recapitulation %",
             "IFNb recapitulation %"),
  Value = c(length(jdms_all),
            length(jdms_up),
            length(jdms_down),
            length(ifna_all),
            length(ifna_up),
            length(ifna_down),
            length(ifnb_all),
            length(ifnb_up),
            length(ifnb_down),
            length(overlap_all_ifna),
            length(overlap_up_ifna),
            length(overlap_down_ifna),
            length(overlap_all_ifnb),
            length(overlap_up_ifnb),
            length(overlap_down_ifnb),
            round(recap_all_ifna, 1),
            round(recap_all_ifnb, 1))
)

print(overlap_summary)
write.csv(overlap_summary, "Overlap_Summary_PatientMuscle_vs_Myobundles.csv", row.names = FALSE)

# Save overlap gene lists (WITH GENE SYMBOLS)
write.table(overlap_all_ifna, "Overlap_JDMS_IFNa_ALL_symbols.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(overlap_up_ifna, "Overlap_JDMS_IFNa_UP_symbols.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(overlap_down_ifna, "Overlap_JDMS_IFNa_DOWN_symbols.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(overlap_all_ifnb, "Overlap_JDMS_IFNb_ALL_symbols.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(overlap_up_ifnb, "Overlap_JDMS_IFNb_UP_symbols.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(overlap_down_ifnb, "Overlap_JDMS_IFNb_DOWN_symbols.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("\nGene lists saved (with gene symbols)\n\n")

# Also save full lists (all significant genes from each comparison)
write.table(jdms_all, "SignificantGenes_JDMS_ALL_symbols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(jdms_up, "SignificantGenes_JDMS_UP_symbols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(jdms_down, "SignificantGenes_JDMS_DOWN_symbols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(ifna_all, "SignificantGenes_IFNa_ALL_symbols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ifna_up, "SignificantGenes_IFNa_UP_symbols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ifna_down, "SignificantGenes_IFNa_DOWN_symbols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(ifnb_all, "SignificantGenes_IFNb_ALL_symbols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ifnb_up, "SignificantGenes_IFNb_UP_symbols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ifnb_down, "SignificantGenes_IFNb_DOWN_symbols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("All significant gene lists saved\n\n")

# Venn diagrams
cat("Creating Venn diagrams...\n")

# 3-way Venn - Upregulated
venn.plot.up <- venn.diagram(
  x = list(JDMS = jdms_up, IFNa = ifna_up, IFNb = ifnb_up),
  filename = NULL,
  output = TRUE,
  main = "Upregulated Genes: JDMS vs IFNa vs IFNb",
  main.cex = 1.5,
  lwd = 2,
  col = c("#C62828", "#1976D2", "#F57C00"),
  fill = c(alpha("#C62828", 0.3), alpha("#1976D2", 0.3), alpha("#F57C00", 0.3)),
  cex = 1.5,
  fontface = "bold",
  cat.cex = 1.3,
  cat.fontface = "bold"
)

pdf("Venn_JDMS_vs_IFNa_vs_IFNb_UP.pdf", width = 10, height = 10)
grid.draw(venn.plot.up)
dev.off()

# 3-way Venn - Downregulated
venn.plot.down <- venn.diagram(
  x = list(JDMS = jdms_down, IFNa = ifna_down, IFNb = ifnb_down),
  filename = NULL,
  output = TRUE,
  main = "Downregulated Genes: JDMS vs IFNa vs IFNb",
  main.cex = 1.5,
  lwd = 2,
  col = c("#C62828", "#1976D2", "#F57C00"),
  fill = c(alpha("#C62828", 0.3), alpha("#1976D2", 0.3), alpha("#F57C00", 0.3)),
  cex = 1.5,
  fontface = "bold",
  cat.cex = 1.3,
  cat.fontface = "bold"
)

pdf("Venn_JDMS_vs_IFNa_vs_IFNb_DOWN.pdf", width = 10, height = 10)
grid.draw(venn.plot.down)
dev.off()

# 3-way Venn - All significant
venn.plot.all <- venn.diagram(
  x = list(JDMS = jdms_all, IFNa = ifna_all, IFNb = ifnb_all),
  filename = NULL,
  output = TRUE,
  main = "All Significant Genes: JDMS vs IFNa vs IFNb",
  main.cex = 1.5,
  lwd = 2,
  col = c("#C62828", "#1976D2", "#F57C00"),
  fill = c(alpha("#C62828", 0.3), alpha("#1976D2", 0.3), alpha("#F57C00", 0.3)),
  cex = 1.5,
  fontface = "bold",
  cat.cex = 1.3,
  cat.fontface = "bold"
)

pdf("Venn_JDMS_vs_IFNa_vs_IFNb_ALL.pdf", width = 10, height = 10)
grid.draw(venn.plot.all)
dev.off()

cat("Venn diagrams saved\n\n")

################################################################################
# ANALYSIS COMPLETE
################################################################################

cat("============================================================\n")
cat("TISSUE-STRATIFIED ANALYSIS COMPLETE!\n")
cat("============================================================\n\n")

cat("Generated files:\n")
cat("PART 1 - Patient Muscle:\n")
cat("  - DESeq2_PatientMuscle_JDMS_vs_Healthy.csv (WITH GENE SYMBOLS)\n")
cat("  - PCA_PatientMuscle_Combined.pdf\n")
cat("  - Heatmap_PatientMuscle_SampleDistances.pdf\n")
cat("  - Volcano_PatientMuscle_JDMS_vs_Healthy.pdf\n\n")

cat("PART 2 - Myobundles:\n")
cat("  - DESeq2_Myobundles_IFNa_vs_NT.csv (WITH GENE SYMBOLS)\n")
cat("  - DESeq2_Myobundles_IFNb_vs_NT.csv (WITH GENE SYMBOLS)\n")
cat("  - DESeq2_Myobundles_IFNb_vs_IFNa.csv (WITH GENE SYMBOLS)\n")
cat("  - PCA_Myobundles_AllDoses.pdf\n")
cat("  - Heatmap_Myobundles_SampleDistances.pdf\n")
cat("  - Volcano_Myobundles_IFNa_vs_NT.pdf\n")
cat("  - Volcano_Myobundles_IFNb_vs_NT.pdf\n\n")

cat("PART 3 - Overlap Analysis:\n")
cat("  - DESeq2_Summary_Statistics.csv\n")
cat("  - Overlap_Summary_PatientMuscle_vs_Myobundles.csv\n")
cat("  - 6 overlap gene lists (TXT files WITH GENE SYMBOLS)\n")
cat("  - 9 complete significant gene lists (TXT files WITH GENE SYMBOLS)\n")
cat("  - 3 Venn diagrams (up, down, all)\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("- IFNa recapitulates %.1f%% of JDMS signature (%d/%d genes)\n", 
            recap_all_ifna, length(overlap_all_ifna), length(jdms_all)))
cat(sprintf("- IFNb recapitulates %.1f%% of JDMS signature (%d/%d genes)\n", 
            recap_all_ifnb, length(overlap_all_ifnb), length(jdms_all)))
cat(sprintf("- Patient muscle: %d significant genes (JDMS vs Healthy)\n", sig_muscle_total))
cat(sprintf("- Myobundles: %d IFNa-responsive genes, %d IFNb-responsive genes\n", 
            sig_ifna_total, sig_ifnb_total))

cat("\n=== ALL RESULTS NOW INCLUDE GENE SYMBOLS ===\n")
