################################################################################
# Extract 10 Specific Contrasts from Combined DESeq2 Object
# Design: ~ Group (with lumped IFN groups)
################################################################################

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

################################################################################
# CHECK THAT REQUIRED OBJECTS EXIST
################################################################################

if(!exists("dds")) {
  stop("ERROR: 'dds' object not found. Please create and run DESeq first.")
}

if(!exists("sample_sheet")) {
  stop("ERROR: 'sample_sheet' not found. Please load your sample metadata.")
}

cat("=== EXTRACTING 10 SPECIFIC CONTRASTS ===\n\n")

# Check what groups we have
cat("Available groups in dds:\n")
print(levels(dds$Group))
cat("\n")

# Check reference level
cat("Reference level:", ref(dds$Group), "\n\n")

################################################################################
# FUNCTION: ADD GENE SYMBOLS
################################################################################

add_gene_symbols <- function(results_df, gene_col = "gene_id") {
  cat("  Adding gene symbols...\n")
  
  ensembl_ids <- results_df[[gene_col]]
  
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  
  results_df$gene_symbol <- ifelse(is.na(gene_symbols), 
                                   results_df[[gene_col]], 
                                   gene_symbols)
  
  col_order <- c(gene_col, "gene_symbol", 
                 setdiff(names(results_df), c(gene_col, "gene_symbol")))
  results_df <- results_df[, col_order]
  
  return(results_df)
}

################################################################################
# DEFINE THE 10 SPECIFIC CONTRASTS
################################################################################

cat("Defining the 10 contrasts...\n\n")

all_contrasts <- list(
  # 1. Disease signature
  "JDM_vs_Healthy" = c("Group", "Muscle-JDM", "Muscle-Healthy"),
  
  # 2-4. Each group vs Healthy muscle
  "Control_vs_Healthy" = c("Group", "Myobundle-Control", "Muscle-Healthy"),
  "INFa_vs_Healthy" = c("Group", "Myobundle-INFa", "Muscle-Healthy"),
  "INFb_vs_Healthy" = c("Group", "Myobundle-INFb", "Muscle-Healthy"),
  
  # 5-7. Each group vs JDM muscle (KEY VALIDATION!)
  "Control_vs_JDM" = c("Group", "Myobundle-Control", "Muscle-JDM"),
  "INFa_vs_JDM" = c("Group", "Myobundle-INFa", "Muscle-JDM"),
  "INFb_vs_JDM" = c("Group", "Myobundle-INFb", "Muscle-JDM"),
  
  # 8-9. IFN treatments vs Control
  "INFa_vs_Control" = c("Group", "Myobundle-INFa", "Myobundle-Control"),
  "INFb_vs_Control" = c("Group", "Myobundle-INFb", "Myobundle-Control"),
  
  # 10. IFNb vs IFNa
  "INFb_vs_INFa" = c("Group", "Myobundle-INFb", "Myobundle-INFa")
)

cat("Total contrasts to extract:", length(all_contrasts), "\n\n")

################################################################################
# EXTRACT ALL RESULTS
################################################################################

cat("Extracting results for all contrasts...\n\n")

results_list <- list()

for(contrast_name in names(all_contrasts)) {
  cat("Extracting:", contrast_name, "...\n")
  
  # Extract results
  res <- results(dds, 
                 contrast = all_contrasts[[contrast_name]],
                 alpha = 0.05)
  
  # Convert to data frame
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", 
                       "lfcSE", "stat", "pvalue", "padj")]
  
  # Add gene symbols
  res_df <- add_gene_symbols(res_df)
  
  # Store in list
  results_list[[contrast_name]] <- res_df
  
  # Print summary
  sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
  sig_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
  sig_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
  
  cat("  Significant:", sig, "(", sig_up, "up,", sig_down, "down)\n\n")
}

################################################################################
# SAVE ALL RESULTS
################################################################################

cat("Saving all results to CSV files...\n")

# Create output directory
dir.create("TenContrasts_Results", showWarnings = FALSE)

for(contrast_name in names(results_list)) {
  filename <- paste0("TenContrasts_Results/DESeq2_", contrast_name, ".csv")
  write.csv(results_list[[contrast_name]], filename, row.names = FALSE)
}

cat("All results saved to TenContrasts_Results/ directory\n\n")

################################################################################
# CREATE SUMMARY TABLE
################################################################################

cat("Creating summary statistics table...\n")

summary_stats <- data.frame(
  Comparison = names(results_list),
  Total_Significant = sapply(results_list, function(x) {
    sum(x$padj < 0.05, na.rm = TRUE)
  }),
  Upregulated = sapply(results_list, function(x) {
    sum(x$padj < 0.05 & x$log2FoldChange > 0, na.rm = TRUE)
  }),
  Downregulated = sapply(results_list, function(x) {
    sum(x$padj < 0.05 & x$log2FoldChange < 0, na.rm = TRUE)
  })
)

print(summary_stats)

write.csv(summary_stats, 
          "TenContrasts_Results/Summary_Statistics.csv", 
          row.names = FALSE)

################################################################################
# EXTRACT SIGNIFICANT GENE LISTS
################################################################################

cat("\nExtracting significant gene lists...\n")

# Function to get significant genes by direction
get_sig_genes <- function(res_df, direction = "both") {
  if(direction == "up") {
    genes <- res_df$gene_symbol[res_df$padj < 0.05 & 
                                 res_df$log2FoldChange > 0 &
                                 !is.na(res_df$padj)]
  } else if(direction == "down") {
    genes <- res_df$gene_symbol[res_df$padj < 0.05 & 
                                 res_df$log2FoldChange < 0 &
                                 !is.na(res_df$padj)]
  } else {
    genes <- res_df$gene_symbol[res_df$padj < 0.05 &
                                 !is.na(res_df$padj)]
  }
  return(genes)
}

# Extract gene lists
gene_lists <- list()

for(contrast_name in names(results_list)) {
  gene_lists[[paste0(contrast_name, "_ALL")]] <- 
    get_sig_genes(results_list[[contrast_name]], "both")
  gene_lists[[paste0(contrast_name, "_UP")]] <- 
    get_sig_genes(results_list[[contrast_name]], "up")
  gene_lists[[paste0(contrast_name, "_DOWN")]] <- 
    get_sig_genes(results_list[[contrast_name]], "down")
}

# Save gene lists
cat("Saving significant gene lists...\n")

for(list_name in names(gene_lists)) {
  if(length(gene_lists[[list_name]]) > 0) {
    filename <- paste0("TenContrasts_Results/GeneList_", list_name, ".txt")
    write.table(gene_lists[[list_name]], 
                filename,
                row.names = FALSE, 
                col.names = FALSE, 
                quote = FALSE)
  }
}

cat("Gene lists saved\n\n")

################################################################################
# CREATE VOLCANO PLOTS FOR ALL 10 COMPARISONS
################################################################################

cat("Creating volcano plots...\n")

# Function to create volcano plot
make_volcano <- function(res_df, title, filename) {
  res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                               "Significant", "Not Significant")
  res_df$significant[is.na(res_df$significant)] <- "Not Significant"
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 P-value") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(filename, p, width = 8, height = 6)
  return(p)
}

# Create volcano plots for all comparisons
for(comp in names(results_list)) {
  make_volcano(results_list[[comp]], 
               comp,
               paste0("TenContrasts_Results/Volcano_", comp, ".pdf"))
}

cat("Volcano plots saved\n\n")

################################################################################
# KEY OVERLAP ANALYSES: JDM RECAPITULATION
################################################################################

cat("=== JDM RECAPITULATION ANALYSIS ===\n\n")

# Get JDM signature
jdm_genes <- get_sig_genes(results_list$JDM_vs_Healthy, "both")
jdm_up <- get_sig_genes(results_list$JDM_vs_Healthy, "up")
jdm_down <- get_sig_genes(results_list$JDM_vs_Healthy, "down")

cat("JDM SIGNATURE:\n")
cat("  Total:", length(jdm_genes), "genes\n")
cat("  Upregulated:", length(jdm_up), "genes\n")
cat("  Downregulated:", length(jdm_down), "genes\n\n")

# Compare each treatment vs Control
treatment_names <- c("INFa_vs_Control", "INFb_vs_Control")

cat("TREATMENT EFFECTS (vs Control):\n\n")

for(treatment in treatment_names) {
  treat_genes <- get_sig_genes(results_list[[treatment]], "both")
  treat_up <- get_sig_genes(results_list[[treatment]], "up")
  treat_down <- get_sig_genes(results_list[[treatment]], "down")
  
  # Overlap with JDM signature
  overlap_total <- length(intersect(jdm_genes, treat_genes))
  overlap_up <- length(intersect(jdm_up, treat_up))
  overlap_down <- length(intersect(jdm_down, treat_down))
  
  pct_total <- (overlap_total / length(jdm_genes)) * 100
  pct_up <- (overlap_up / length(jdm_up)) * 100
  pct_down <- (overlap_down / length(jdm_down)) * 100
  
  cat(treatment, "recapitulates JDM:\n")
  cat(sprintf("  Total: %d/%d genes (%.1f%%)\n", 
              overlap_total, length(jdm_genes), pct_total))
  cat(sprintf("  Up: %d/%d genes (%.1f%%)\n", 
              overlap_up, length(jdm_up), pct_up))
  cat(sprintf("  Down: %d/%d genes (%.1f%%)\n\n", 
              overlap_down, length(jdm_down), pct_down))
}

################################################################################
# CREATE VENN DIAGRAMS
################################################################################

cat("Creating Venn diagrams...\n")

# Function to create Venn diagram
make_venn <- function(set1, set2, name1, name2, title, filename) {
  if(length(set1) == 0 || length(set2) == 0) {
    cat("  Skipping", title, "- one or both sets empty\n")
    return(NULL)
  }
  
  venn.plot <- venn.diagram(
    x = list(set1, set2),
    category.names = c(name1, name2),
    filename = NULL,
    output = TRUE,
    main = title,
    main.cex = 1.5,
    lwd = 2,
    col = c("#C62828", "#1976D2"),
    fill = c(alpha("#C62828", 0.3), alpha("#1976D2", 0.3)),
    cex = 1.5,
    fontface = "bold",
    cat.cex = 1.3,
    cat.fontface = "bold"
  )
  
  pdf(filename, width = 8, height = 8)
  grid.draw(venn.plot)
  dev.off()
}

# JDM vs INFa
make_venn(jdm_up, get_sig_genes(results_list$INFa_vs_Control, "up"),
          "JDM", "IFNa",
          "Upregulated Genes: JDM vs IFNa",
          "TenContrasts_Results/Venn_JDM_vs_INFa_UP.pdf")

make_venn(jdm_down, get_sig_genes(results_list$INFa_vs_Control, "down"),
          "JDM", "IFNa",
          "Downregulated Genes: JDM vs IFNa",
          "TenContrasts_Results/Venn_JDM_vs_INFa_DOWN.pdf")

# JDM vs INFb
make_venn(jdm_up, get_sig_genes(results_list$INFb_vs_Control, "up"),
          "JDM", "IFNb",
          "Upregulated Genes: JDM vs IFNb",
          "TenContrasts_Results/Venn_JDM_vs_INFb_UP.pdf")

make_venn(jdm_down, get_sig_genes(results_list$INFb_vs_Control, "down"),
          "JDM", "IFNb",
          "Downregulated Genes: JDM vs IFNb",
          "TenContrasts_Results/Venn_JDM_vs_INFb_DOWN.pdf")

# 3-way Venn: JDM vs INFa vs INFb (upregulated)
if(length(jdm_up) > 0 && 
   length(get_sig_genes(results_list$INFa_vs_Control, "up")) > 0 &&
   length(get_sig_genes(results_list$INFb_vs_Control, "up")) > 0) {
  
  venn.plot.3way <- venn.diagram(
    x = list(JDM = jdm_up, 
             INFa = get_sig_genes(results_list$INFa_vs_Control, "up"),
             INFb = get_sig_genes(results_list$INFb_vs_Control, "up")),
    filename = NULL,
    output = TRUE,
    main = "Upregulated Genes: JDM vs INFa vs INFb",
    main.cex = 1.5,
    lwd = 2,
    col = c("#C62828", "#1976D2", "#F57C00"),
    fill = c(alpha("#C62828", 0.3), alpha("#1976D2", 0.3), alpha("#F57C00", 0.3)),
    cex = 1.5,
    fontface = "bold",
    cat.cex = 1.3,
    cat.fontface = "bold"
  )
  
  pdf("TenContrasts_Results/Venn_JDM_vs_INFa_vs_INFb_UP.pdf", width = 10, height = 10)
  grid.draw(venn.plot.3way)
  dev.off()
}

cat("Venn diagrams saved\n\n")

################################################################################
# CORRELATION ANALYSIS
################################################################################

cat("=== CORRELATION ANALYSIS ===\n\n")

# Merge JDM and INFa results for correlation
merged_ifna <- merge(
  results_list$JDM_vs_Healthy[, c("gene_id", "gene_symbol", "log2FoldChange", "padj")],
  results_list$INFa_vs_Control[, c("gene_id", "log2FoldChange", "padj")],
  by = "gene_id",
  suffixes = c("_JDM", "_INFa")
)

merged_ifna <- merged_ifna[!is.na(merged_ifna$log2FoldChange_JDM) & 
                           !is.na(merged_ifna$log2FoldChange_INFa), ]

# Filter for significant genes
merged_ifna_sig <- merged_ifna[merged_ifna$padj_JDM < 0.05 | merged_ifna$padj_INFa < 0.05, ]

if(nrow(merged_ifna_sig) > 0) {
  cor_ifna <- cor(merged_ifna_sig$log2FoldChange_JDM, merged_ifna_sig$log2FoldChange_INFa)
  cat("Correlation (JDM vs INFa, sig genes):", round(cor_ifna, 3), "\n")
  
  # Scatter plot
  p_ifna <- ggplot(merged_ifna_sig, aes(x = log2FoldChange_JDM, y = log2FoldChange_INFa)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(title = "JDM vs IFNa: Log2 Fold Change Correlation",
         subtitle = paste0("Pearson r = ", round(cor_ifna, 3)),
         x = "Log2 Fold Change (JDM vs Healthy)",
         y = "Log2 Fold Change (IFNa vs Control)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  
  ggsave("TenContrasts_Results/Scatter_JDM_vs_INFa.pdf", p_ifna, width = 8, height = 8)
}

# Same for INFb
merged_ifnb <- merge(
  results_list$JDM_vs_Healthy[, c("gene_id", "gene_symbol", "log2FoldChange", "padj")],
  results_list$INFb_vs_Control[, c("gene_id", "log2FoldChange", "padj")],
  by = "gene_id",
  suffixes = c("_JDM", "_INFb")
)

merged_ifnb <- merged_ifnb[!is.na(merged_ifnb$log2FoldChange_JDM) & 
                           !is.na(merged_ifnb$log2FoldChange_INFb), ]

merged_ifnb_sig <- merged_ifnb[merged_ifnb$padj_JDM < 0.05 | merged_ifnb$padj_INFb < 0.05, ]

if(nrow(merged_ifnb_sig) > 0) {
  cor_ifnb <- cor(merged_ifnb_sig$log2FoldChange_JDM, merged_ifnb_sig$log2FoldChange_INFb)
  cat("Correlation (JDM vs INFb, sig genes):", round(cor_ifnb, 3), "\n\n")
  
  # Scatter plot
  p_ifnb <- ggplot(merged_ifnb_sig, aes(x = log2FoldChange_JDM, y = log2FoldChange_INFb)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(title = "JDM vs IFNb: Log2 Fold Change Correlation",
         subtitle = paste0("Pearson r = ", round(cor_ifnb, 3)),
         x = "Log2 Fold Change (JDM vs Healthy)",
         y = "Log2 Fold Change (IFNb vs Control)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  
  ggsave("TenContrasts_Results/Scatter_JDM_vs_INFb.pdf", p_ifnb, width = 8, height = 8)
}

cat("Correlation plots saved\n\n")

################################################################################
# ANALYSIS COMPLETE
################################################################################

cat("============================================================\n")
cat("10 CONTRASTS EXTRACTED AND ANALYZED!\n")
cat("============================================================\n\n")

cat("Generated files in TenContrasts_Results/:\n")
cat("- 10 CSV files with full DESeq2 results (with gene symbols)\n")
cat("- 30 gene lists (all/up/down for each comparison)\n")
cat("- 10 volcano plots\n")
cat("- 4-5 Venn diagrams\n")
cat("- 2 correlation scatter plots\n")
cat("- Summary_Statistics.csv\n\n")

cat("Objects in workspace:\n")
cat("- results_list: All DESeq2 results\n")
cat("- gene_lists: All significant gene lists\n")
cat("- summary_stats: Summary table\n\n")

cat("KEY COMPARISONS FOR YOUR RESEARCH QUESTION:\n")
cat("1. JDM_vs_Healthy - Disease signature\n")
cat("2. INFa_vs_Control - IFNa treatment effect\n")
cat("3. INFb_vs_Control - IFNb treatment effect\n")
cat("4. INFa_vs_JDM - Does IFNa recapitulate JDM?\n")
cat("5. INFb_vs_JDM - Does IFNb recapitulate JDM?\n\n")

cat("NEXT STEPS:\n")
cat("1. Run GSEA on these 10 comparisons\n")
cat("2. Identify top overlapping genes between JDM and treatments\n")
cat("3. Determine which IFN type is the better model\n")
