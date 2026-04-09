# =============================================================================
# TOP 100 VARIABLE GENES HEATMAP - FIXED FOR EMPTY PDF
# =============================================================================
# Purpose: Show expression patterns of most variable genes across samples
# Input: dds object already in environment
# Output: PDF heatmap with hierarchical clustering
# =============================================================================

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(genefilter)  # <-- NEEDED FOR rowVars()

cat("Starting top 100 variable genes analysis...\n")

# -----------------------------------------------------------------------------
# 1. TRANSFORM DATA
# -----------------------------------------------------------------------------
cat("Performing variance stabilizing transformation...\n")
vsd <- vst(dds, blind = FALSE)

# -----------------------------------------------------------------------------
# 2. SELECT TOP 100 VARIABLE GENES
# -----------------------------------------------------------------------------
cat("Calculating gene variances...\n")

# Calculate variance for each gene
gene_vars <- rowVars(assay(vsd))
topVarGenes <- head(order(gene_vars, decreasing = TRUE), 100)

cat(paste0("Selected top 100 genes (variance range: ", 
           round(min(gene_vars[topVarGenes]), 3), " to ", 
           round(max(gene_vars[topVarGenes]), 3), ")\n"))

# Extract the data for these genes
mat <- assay(vsd)[topVarGenes, ]
cat(paste0("Matrix dimensions: ", nrow(mat), " genes x ", ncol(mat), " samples\n"))

# -----------------------------------------------------------------------------
# 3. SCALE THE DATA (Z-score by row)
# -----------------------------------------------------------------------------
cat("Scaling data by gene (Z-score)...\n")
mat_scaled <- t(scale(t(mat)))

# Check for any issues
cat(paste0("Scaled matrix has ", sum(is.na(mat_scaled)), " NA values\n"))

# -----------------------------------------------------------------------------
# 4. CREATE ANNOTATION FOR SAMPLES
# -----------------------------------------------------------------------------
cat("\nAvailable colData columns:\n")
print(colnames(colData(dds)))

# Build annotation only with columns that exist and have no NAs
annotation_col <- data.frame(row.names = colnames(dds))

# Try to add each column if it exists
if("Tissue" %in% colnames(colData(dds))) {
  annotation_col$Tissue <- as.factor(colData(dds)$Tissue)
}

if("Condition" %in% colnames(colData(dds))) {
  annotation_col$Condition <- as.factor(colData(dds)$Condition)
}

if("Location" %in% colnames(colData(dds))) {
  annotation_col$Location <- as.factor(colData(dds)$Location)
}

# Remove any columns with NAs
annotation_col <- annotation_col[, colSums(is.na(annotation_col)) == 0, drop = FALSE]

cat("\nAnnotation columns used:\n")
print(head(annotation_col))

# -----------------------------------------------------------------------------
# 5. GENERATE HEATMAP
# -----------------------------------------------------------------------------
cat("\nGenerating heatmap...\n")

# Color scheme: blue (low) -> white (medium) -> red (high)
colors <- colorRampPalette(c("blue", "white", "red"))(255)

pdf("top100_variable_genes_heatmap.pdf", width = 14, height = 12)

# If annotation exists, use it. Otherwise, plot without annotation
if(ncol(annotation_col) > 0) {
  pheatmap(
    mat_scaled,
    color = colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,  # Too many genes to show names
    show_colnames = TRUE,
    annotation_col = annotation_col,
    main = "Top 100 Most Variable Genes",
    fontsize = 10,
    fontsize_col = 8,
    scale = "none"  # Already scaled above
  )
} else {
  # Plot without annotation if columns don't exist
  pheatmap(
    mat_scaled,
    color = colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    main = "Top 100 Most Variable Genes",
    fontsize = 10,
    fontsize_col = 8,
    scale = "none"
  )
}

dev.off()

cat("✅ Top 100 variable genes heatmap saved to: top100_variable_genes_heatmap.pdf\n")

# -----------------------------------------------------------------------------
# 6. Save gene list
# -----------------------------------------------------------------------------
top_gene_names <- rownames(mat)
write.csv(
  data.frame(
    Gene = top_gene_names,
    Variance = gene_vars[topVarGenes]
  ),
  "top100_variable_genes_list.csv",
  row.names = FALSE
)

cat("✅ Gene list saved to: top100_variable_genes_list.csv\n")
cat("\n=== DONE ===\n")
