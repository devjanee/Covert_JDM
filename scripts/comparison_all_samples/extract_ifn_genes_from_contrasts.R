################################################################################
# Extract IFN Genes from All 10 Contrasts
# Reads IFN gene list and pulls out their stats from each comparison
################################################################################

library(readxl)
library(dplyr)
library(tidyr)

################################################################################
# 1. LOAD IFN GENES LIST
################################################################################

cat("=== EXTRACTING IFN GENES FROM ALL CONTRASTS ===\n\n")

# Option 1: Read from Excel file
cat("Reading IFN genes from IFN_genes.xlsx...\n")

tryCatch({
  ifn_genes_df <- read_excel("~/Documents/Covert_11358/IFN_genes.xlsx")
  cat("Successfully read", nrow(ifn_genes_df), "genes from Excel file\n\n")
}, error = function(e) {
  cat("Could not read IFN_genes.xlsx\n")
  cat("Error:", e$message, "\n\n")
  cat("Please provide your IFN genes list!\n")
  cat("Either:\n")
  cat("1. Update the 'ifn_gene_list' vector below, or\n")
  cat("2. Make sure IFN_genes.xlsx is in your working directory\n\n")
  stop("IFN genes file not found")
})

# Extract gene symbols (assuming first column contains gene names)
ifn_gene_list <- ifn_genes_df[[1]]

# Remove any NA values
ifn_gene_list <- ifn_gene_list[!is.na(ifn_gene_list)]

cat("IFN genes to extract:\n")
print(head(ifn_gene_list, 20))
cat("... (showing first 20)\n")
cat("Total IFN genes:", length(ifn_gene_list), "\n\n")

################################################################################
# 2. CHECK THAT RESULTS EXIST
################################################################################

if(!exists("results_list")) {
  cat("ERROR: results_list not found!\n")
  cat("Please run extract_10_contrasts.R first\n\n")
  stop("Results not available")
}

cat("Found results for", length(results_list), "contrasts\n\n")

################################################################################
# 3. EXTRACT IFN GENES FROM EACH CONTRAST
################################################################################

cat("Extracting IFN genes from each contrast...\n\n")

# Create list to store IFN gene results
ifn_results_list <- list()

for(contrast_name in names(results_list)) {
  cat("Processing:", contrast_name, "...\n")
  
  # Get results for this contrast
  res_df <- results_list[[contrast_name]]
  
  # Filter for IFN genes only
  ifn_subset <- res_df[res_df$gene_symbol %in% ifn_gene_list, ]
  
  # Add contrast name column
  ifn_subset$contrast <- contrast_name
  
  # Store
  ifn_results_list[[contrast_name]] <- ifn_subset
  
  cat("  Found", nrow(ifn_subset), "IFN genes\n")
}

cat("\n")

################################################################################
# 4. CREATE COMPREHENSIVE IFN GENES TABLE
################################################################################

cat("Creating comprehensive IFN genes table...\n")

# Combine all contrasts
all_ifn_results <- bind_rows(ifn_results_list)

# Reorder columns for clarity
all_ifn_results <- all_ifn_results %>%
  select(contrast, gene_symbol, gene_id, baseMean, log2FoldChange, 
         lfcSE, stat, pvalue, padj)

# Sort by gene symbol, then contrast
all_ifn_results <- all_ifn_results %>%
  arrange(gene_symbol, contrast)

cat("Total rows:", nrow(all_ifn_results), "\n")
cat("Unique genes found:", length(unique(all_ifn_results$gene_symbol)), "/", 
    length(ifn_gene_list), "\n\n")

################################################################################
# 5. CREATE WIDE FORMAT TABLE (EXCEL-LIKE)
################################################################################

cat("Creating wide format table (one row per gene)...\n")

# For each stat, create wide format
create_wide_table <- function(data, value_col, prefix) {
  data %>%
    select(gene_symbol, contrast, all_of(value_col)) %>%
    pivot_wider(names_from = contrast, 
                values_from = all_of(value_col),
                names_prefix = paste0(prefix, "_"))
}

# Create separate tables for each metric
wide_log2fc <- create_wide_table(all_ifn_results, "log2FoldChange", "log2FC")
wide_padj <- create_wide_table(all_ifn_results, "padj", "padj")
wide_baseMean <- create_wide_table(all_ifn_results, "baseMean", "baseMean")

# Merge all together
ifn_wide <- wide_log2fc %>%
  left_join(wide_padj, by = "gene_symbol") %>%
  left_join(wide_baseMean, by = "gene_symbol")

# Reorder columns nicely
# Group by contrast: log2FC, padj, baseMean for each contrast
contrast_names <- names(results_list)
col_order <- c("gene_symbol")

for(contrast in contrast_names) {
  col_order <- c(col_order,
                 paste0("log2FC_", contrast),
                 paste0("padj_", contrast),
                 paste0("baseMean_", contrast))
}

# Keep only columns that exist
col_order <- col_order[col_order %in% names(ifn_wide)]
ifn_wide <- ifn_wide[, col_order]

cat("Wide table created:", nrow(ifn_wide), "genes x", ncol(ifn_wide), "columns\n\n")

################################################################################
# 6. CREATE SUMMARY TABLE
################################################################################

cat("Creating summary table...\n")

# Count significant IFN genes in each contrast
ifn_summary <- data.frame(
  Contrast = names(results_list),
  Total_IFN_Genes = sapply(ifn_results_list, nrow),
  Significant_IFN = sapply(ifn_results_list, function(x) {
    sum(x$padj < 0.05, na.rm = TRUE)
  }),
  Upregulated_IFN = sapply(ifn_results_list, function(x) {
    sum(x$padj < 0.05 & x$log2FoldChange > 0, na.rm = TRUE)
  }),
  Downregulated_IFN = sapply(ifn_results_list, function(x) {
    sum(x$padj < 0.05 & x$log2FoldChange < 0, na.rm = TRUE)
  })
)

print(ifn_summary)
cat("\n")

################################################################################
# 7. SAVE ALL OUTPUTS
################################################################################

cat("Saving results...\n")

# Create output directory
dir.create("IFN_Genes_Results", showWarnings = FALSE)

# 1. Long format (all data)
write.csv(all_ifn_results, 
          "IFN_Genes_Results/IFN_genes_all_contrasts_LONG.csv", 
          row.names = FALSE)

# 2. Wide format (Excel-like, one row per gene)
write.csv(ifn_wide, 
          "IFN_Genes_Results/IFN_genes_all_contrasts_WIDE.csv", 
          row.names = FALSE)

# 3. Summary table
write.csv(ifn_summary, 
          "IFN_Genes_Results/IFN_genes_summary.csv", 
          row.names = FALSE)

# 4. Separate CSV for each contrast (like original format)
for(contrast_name in names(ifn_results_list)) {
  contrast_data <- ifn_results_list[[contrast_name]] %>%
    select(gene_symbol, gene_id, baseMean, log2FoldChange, 
           lfcSE, stat, pvalue, padj)
  
  write.csv(contrast_data, 
            paste0("IFN_Genes_Results/IFN_genes_", contrast_name, ".csv"),
            row.names = FALSE)
}

cat("\nAll files saved to IFN_Genes_Results/\n\n")

################################################################################
# 8. CREATE HEATMAP OF IFN GENE LOG2FC
################################################################################

cat("Creating heatmap of IFN gene responses...\n")

library(pheatmap)

# Create matrix for heatmap (log2FC only)
heatmap_data <- all_ifn_results %>%
  select(gene_symbol, contrast, log2FoldChange) %>%
  pivot_wider(names_from = contrast, values_from = log2FoldChange) %>%
  as.data.frame()

rownames(heatmap_data) <- heatmap_data$gene_symbol
heatmap_data$gene_symbol <- NULL

# Remove genes with all NA
heatmap_data <- heatmap_data[rowSums(is.na(heatmap_data)) < ncol(heatmap_data), ]

if(nrow(heatmap_data) > 0) {
  # Create heatmap
  pdf("IFN_Genes_Results/IFN_genes_heatmap.pdf", width = 12, height = 10)
  pheatmap(as.matrix(heatmap_data),
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-5, 5, length.out = 101),
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           main = "IFN Genes: Log2 Fold Change Across All Contrasts",
           fontsize_row = 7,
           fontsize_col = 9,
           na_col = "grey90")
  dev.off()
  
  cat("Heatmap saved\n")
} else {
  cat("Not enough data for heatmap\n")
}

################################################################################
# 9. IDENTIFY TOP RESPONSIVE IFN GENES
################################################################################

cat("\n=== TOP RESPONSIVE IFN GENES ===\n\n")

# Find genes that are significant in multiple contrasts
gene_sig_counts <- all_ifn_results %>%
  filter(padj < 0.05) %>%
  group_by(gene_symbol) %>%
  summarize(
    n_significant = n(),
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    max_log2FC = max(abs(log2FoldChange), na.rm = TRUE),
    contrasts = paste(contrast, collapse = "; ")
  ) %>%
  arrange(desc(n_significant), desc(max_log2FC))

cat("IFN genes significant in multiple contrasts:\n")
print(head(gene_sig_counts, 20))

write.csv(gene_sig_counts, 
          "IFN_Genes_Results/IFN_genes_top_responsive.csv",
          row.names = FALSE)

################################################################################
# COMPLETE
################################################################################

cat("\n============================================================\n")
cat("IFN GENE EXTRACTION COMPLETE!\n")
cat("============================================================\n\n")

cat("Generated files in IFN_Genes_Results/:\n")
cat("- IFN_genes_all_contrasts_LONG.csv (all data, long format)\n")
cat("- IFN_genes_all_contrasts_WIDE.csv (Excel-like, one row per gene)\n")
cat("- IFN_genes_summary.csv (counts per contrast)\n")
cat("- IFN_genes_[contrast].csv (10 separate files, one per contrast)\n")
cat("- IFN_genes_heatmap.pdf (visualization)\n")
cat("- IFN_genes_top_responsive.csv (genes significant in multiple contrasts)\n\n")

cat("Key findings:\n")
cat("- Total IFN genes analyzed:", length(ifn_gene_list), "\n")
cat("- IFN genes found in results:", length(unique(all_ifn_results$gene_symbol)), "\n")
cat("- Genes significant in at least one contrast:", 
    nrow(gene_sig_counts), "\n\n")

cat("Most responsive IFN gene:", 
    as.character(gene_sig_counts$gene_symbol[1]), 
    "(significant in", gene_sig_counts$n_significant[1], "contrasts)\n")
