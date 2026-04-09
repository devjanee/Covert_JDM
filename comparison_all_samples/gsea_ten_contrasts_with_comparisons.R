################################################################################
# GSEA Analysis for 10 Contrasts
# Purpose: Run GSEA on all comparisons and create comparison plots
# Style: Similar to tissue_stratified_gsea_with_comparisons.R
################################################################################

library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)

################################################################################
# STEP 1: LOAD HALLMARK GENE SETS
################################################################################

cat("=== GSEA ANALYSIS FOR 10 CONTRASTS ===\n\n")

cat("Loading Hallmark gene sets from MSigDB...\n")
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- hallmark_sets %>% 
  split(x = .$gene_symbol, f = .$gs_name)

cat("Loaded", length(hallmark_list), "Hallmark gene sets\n\n")

################################################################################
# STEP 2: CHECK THAT RESULTS EXIST
################################################################################

if(!exists("results_list")) {
  stop("ERROR: results_list not found! Please run extract_10_contrasts.R first")
}

cat("Found results for", length(results_list), "contrasts\n")
cat("Contrasts:\n")
for(i in seq_along(names(results_list))) {
  cat(sprintf("  %d. %s\n", i, names(results_list)[i]))
}
cat("\n")

################################################################################
# STEP 3: PREPARE RANKED GENE LISTS
################################################################################

cat("Preparing ranked gene lists for GSEA...\n")

# Function to create ranked list from results with gene symbols
prepare_ranked_list <- function(results_df) {
  # Remove rows with NA values
  results_df <- results_df[!is.na(results_df$pvalue) & 
                           !is.na(results_df$log2FoldChange) &
                           !is.na(results_df$gene_symbol), ]
  
  # Cap very small p-values to prevent Inf values
  # Minimum p-value = 1e-300 (corresponds to -log10(p) = 300)
  results_df$pvalue <- pmax(results_df$pvalue, 1e-300)
  
  # Calculate ranking metric: sign(log2FC) * -log10(pvalue)
  results_df$rank_metric <- sign(results_df$log2FoldChange) * -log10(results_df$pvalue)
  
  # Remove any remaining infinite or NA values
  results_df <- results_df[is.finite(results_df$rank_metric), ]
  
  # Handle duplicates - keep the one with best rank metric
  results_df <- results_df %>%
    group_by(gene_symbol) %>%
    slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    as.data.frame()
  
  # Create named vector
  ranks <- results_df$rank_metric
  names(ranks) <- results_df$gene_symbol
  
  # Sort by rank (descending)
  ranks <- sort(ranks, decreasing = TRUE)
  
  return(ranks)
}

# Create ranked lists for all contrasts
ranked_lists <- list()

for(contrast_name in names(results_list)) {
  ranked_lists[[contrast_name]] <- prepare_ranked_list(results_list[[contrast_name]])
  cat(contrast_name, ":", length(ranked_lists[[contrast_name]]), "genes\n")
}
cat("\n")

################################################################################
# STEP 4: RUN FGSEA ON ALL CONTRASTS
################################################################################

cat("Running GSEA on all 10 comparisons...\n\n")

gsea_results <- list()

for(comparison in names(ranked_lists)) {
  cat("Running GSEA:", comparison, "...\n")
  
  # Run fgsea
  fgsea_res <- fgsea(
    pathways = hallmark_list,
    stats = ranked_lists[[comparison]],
    minSize = 15,
    maxSize = 500
  )
  
  # Sort by NES
  fgsea_res <- fgsea_res[order(fgsea_res$NES, decreasing = TRUE), ]
  
  # Store results
  gsea_results[[comparison]] <- fgsea_res
  
  # Print summary
  sig_up <- sum(fgsea_res$padj < 0.05 & fgsea_res$NES > 0, na.rm = TRUE)
  sig_down <- sum(fgsea_res$padj < 0.05 & fgsea_res$NES < 0, na.rm = TRUE)
  cat("   Significant pathways:", sig_up + sig_down, 
      "(", sig_up, "up,", sig_down, "down)\n")
}
cat("\n")

################################################################################
# STEP 5: SAVE GSEA RESULTS
################################################################################

cat("Saving GSEA results to CSV files...\n")

# Create output directory
dir.create("TenContrasts_GSEA", showWarnings = FALSE)

for(comparison in names(gsea_results)) {
  res_df <- as.data.frame(gsea_results[[comparison]])
  res_df$leadingEdge <- sapply(res_df$leadingEdge, 
                               function(x) paste(x, collapse = ";"))
  write.csv(res_df, 
            file = paste0("TenContrasts_GSEA/GSEA_", comparison, ".csv"),
            row.names = FALSE)
}

cat("GSEA results saved\n\n")

################################################################################
# STEP 6: PATHWAY OVERLAP ANALYSIS
################################################################################

cat("=== PATHWAY OVERLAP ANALYSIS ===\n\n")

# Get significant pathways from each comparison
sig_pathways_list <- list()
sig_up_pathways_list <- list()
sig_down_pathways_list <- list()

for(comp in names(gsea_results)) {
  sig_pathways_list[[comp]] <- gsea_results[[comp]]$pathway[
    gsea_results[[comp]]$padj < 0.05
  ]
  sig_up_pathways_list[[comp]] <- gsea_results[[comp]]$pathway[
    gsea_results[[comp]]$padj < 0.05 & gsea_results[[comp]]$NES > 0
  ]
  sig_down_pathways_list[[comp]] <- gsea_results[[comp]]$pathway[
    gsea_results[[comp]]$padj < 0.05 & gsea_results[[comp]]$NES < 0
  ]
}

# Print summary
cat("Significant pathways per contrast:\n")
for(comp in names(sig_pathways_list)) {
  cat(sprintf("  %-20s: %3d total (%2d up, %2d down)\n", 
              comp, 
              length(sig_pathways_list[[comp]]),
              length(sig_up_pathways_list[[comp]]),
              length(sig_down_pathways_list[[comp]])))
}
cat("\n")

# KEY COMPARISON: JDM vs IFN treatments
cat("=== KEY COMPARISON: JDM PATHWAY RECAPITULATION ===\n\n")

jdm_pathways <- sig_pathways_list$JDM_vs_Healthy
ifna_pathways <- sig_pathways_list$INFa_vs_Control
ifnb_pathways <- sig_pathways_list$INFb_vs_Control

overlap_jdm_ifna <- intersect(jdm_pathways, ifna_pathways)
overlap_jdm_ifnb <- intersect(jdm_pathways, ifnb_pathways)

cat("JDM signature pathways:", length(jdm_pathways), "\n")
cat("INFa treatment pathways:", length(ifna_pathways), "\n")
cat("INFb treatment pathways:", length(ifnb_pathways), "\n\n")

cat("JDM-INFa pathway overlap:", length(overlap_jdm_ifna), 
    sprintf("(%.1f%% of JDM)\n", 100 * length(overlap_jdm_ifna) / max(length(jdm_pathways), 1)))
cat("JDM-INFb pathway overlap:", length(overlap_jdm_ifnb), 
    sprintf("(%.1f%% of JDM)\n\n", 100 * length(overlap_jdm_ifnb) / max(length(jdm_pathways), 1)))

if(length(overlap_jdm_ifna) > 0) {
  cat("CONSENSUS PATHWAYS (JDM-INFa):\n")
  for(p in overlap_jdm_ifna) {
    clean_name <- gsub("HALLMARK_", "", p)
    clean_name <- gsub("_", " ", clean_name)
    cat("  -", clean_name, "\n")
  }
  cat("\n")
}

################################################################################
# STEP 7: VENN DIAGRAMS
################################################################################

cat("Creating pathway overlap Venn diagrams...\n")

# 3-way Venn: JDM vs INFa vs INFb (all significant pathways)
if(length(jdm_pathways) > 0 && length(ifna_pathways) > 0 && length(ifnb_pathways) > 0) {
  venn.plot.all <- venn.diagram(
    x = list(JDM = jdm_pathways, 
             INFa = ifna_pathways, 
             INFb = ifnb_pathways),
    category.names = c("JDM\n(Muscle)", "IFNa\n(Treatment)", "IFNb\n(Treatment)"),
    filename = NULL,
    output = TRUE,
    main = "Pathway Overlap: JDM vs IFN Treatments",
    main.cex = 1.5,
    lwd = 2,
    col = c("#C62828", "#1976D2", "#F57C00"),
    fill = c(alpha("#C62828", 0.3), alpha("#1976D2", 0.3), alpha("#F57C00", 0.3)),
    cex = 1.5,
    fontface = "bold",
    cat.cex = 1.3,
    cat.fontface = "bold"
  )
  
  pdf("TenContrasts_GSEA/Venn_JDM_vs_INFa_vs_INFb_Pathways.pdf", width = 10, height = 10)
  grid.draw(venn.plot.all)
  dev.off()
}

# Upregulated pathways
jdm_up <- sig_up_pathways_list$JDM_vs_Healthy
ifna_up <- sig_up_pathways_list$INFa_vs_Control
ifnb_up <- sig_up_pathways_list$INFb_vs_Control

if(length(jdm_up) > 0 && length(ifna_up) > 0 && length(ifnb_up) > 0) {
  venn.plot.up <- venn.diagram(
    x = list(JDM = jdm_up, 
             INFa = ifna_up, 
             INFb = ifnb_up),
    category.names = c("JDM\n(Muscle)", "IFNa\n(Treatment)", "IFNb\n(Treatment)"),
    filename = NULL,
    output = TRUE,
    main = "Upregulated Pathways: JDM vs IFN Treatments",
    main.cex = 1.5,
    lwd = 2,
    col = c("#C62828", "#1976D2", "#F57C00"),
    fill = c(alpha("#C62828", 0.3), alpha("#1976D2", 0.3), alpha("#F57C00", 0.3)),
    cex = 1.5,
    fontface = "bold",
    cat.cex = 1.3,
    cat.fontface = "bold"
  )
  
  pdf("TenContrasts_GSEA/Venn_JDM_vs_INFa_vs_INFb_Pathways_UP.pdf", width = 10, height = 10)
  grid.draw(venn.plot.up)
  dev.off()
}

cat("Venn diagrams saved\n\n")

################################################################################
# STEP 8: CORRELATION OF NES SCORES
################################################################################

cat("Analyzing correlation of pathway enrichment scores...\n")

# Merge JDM, INFa, and INFb GSEA results
merged_gsea <- merge(
  gsea_results$JDM_vs_Healthy[, c("pathway", "NES", "padj")],
  gsea_results$INFa_vs_Control[, c("pathway", "NES", "padj")],
  by = "pathway",
  suffixes = c("_JDM", "_INFa")
)

merged_gsea <- merge(
  merged_gsea,
  gsea_results$INFb_vs_Control[, c("pathway", "NES", "padj")],
  by = "pathway"
)
names(merged_gsea)[names(merged_gsea) == "NES"] <- "NES_INFb"
names(merged_gsea)[names(merged_gsea) == "padj"] <- "padj_INFb"

# Calculate correlations
cor_jdm_ifna <- cor(merged_gsea$NES_JDM, merged_gsea$NES_INFa)
cor_jdm_ifnb <- cor(merged_gsea$NES_JDM, merged_gsea$NES_INFb)
cor_ifna_ifnb <- cor(merged_gsea$NES_INFa, merged_gsea$NES_INFb)

cat("Correlation of NES (all pathways):\n")
cat("  JDM vs INFa:", round(cor_jdm_ifna, 3), "\n")
cat("  JDM vs INFb:", round(cor_jdm_ifnb, 3), "\n")
cat("  INFa vs INFb:", round(cor_ifna_ifnb, 3), "\n\n")

# Scatter plots
p_jdm_ifna <- ggplot(merged_gsea, aes(x = NES_JDM, y = NES_INFa)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Pathway Enrichment: JDM vs IFNa",
       subtitle = paste0("Pearson r = ", round(cor_jdm_ifna, 3)),
       x = "NES (JDM - Muscle Disease)",
       y = "NES (IFNa - Treatment)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("TenContrasts_GSEA/Scatter_JDM_vs_INFa_Pathway_NES.pdf", 
       p_jdm_ifna, width = 8, height = 8)

p_jdm_ifnb <- ggplot(merged_gsea, aes(x = NES_JDM, y = NES_INFb)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Pathway Enrichment: JDM vs IFNb",
       subtitle = paste0("Pearson r = ", round(cor_jdm_ifnb, 3)),
       x = "NES (JDM - Muscle Disease)",
       y = "NES (IFNb - Treatment)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("TenContrasts_GSEA/Scatter_JDM_vs_INFb_Pathway_NES.pdf", 
       p_jdm_ifnb, width = 8, height = 8)

cat("Correlation plots saved\n\n")

################################################################################
# STEP 9: SIDE-BY-SIDE COMPARISON BARPLOTS
################################################################################

cat("Creating side-by-side pathway comparison barplots...\n")

# Function to create comparison barplot
create_comparison_barplot <- function(gsea1, gsea2, name1, name2, 
                                     top_n = 15, filename = NULL) {
  # Get top pathways from first comparison (by absolute NES)
  top_pathways <- head(gsea1[order(-abs(gsea1$NES)), ], top_n)$pathway
  
  # Get NES for these pathways in both comparisons
  comparison_df <- data.frame(
    pathway = top_pathways,
    NES1 = gsea1$NES[match(top_pathways, gsea1$pathway)],
    padj1 = gsea1$padj[match(top_pathways, gsea1$pathway)],
    NES2 = gsea2$NES[match(top_pathways, gsea2$pathway)],
    padj2 = gsea2$padj[match(top_pathways, gsea2$pathway)]
  )
  
  # Clean pathway names
  comparison_df$pathway_clean <- gsub("HALLMARK_", "", comparison_df$pathway)
  comparison_df$pathway_clean <- gsub("_", " ", comparison_df$pathway_clean)
  
  # Reshape for plotting
  comparison_long <- comparison_df[, c("pathway_clean", "NES1", "NES2")]
  names(comparison_long)[2:3] <- c(name1, name2)
  
  comparison_long <- pivot_longer(comparison_long,
                                  cols = c(all_of(name1), all_of(name2)),
                                  names_to = "Comparison",
                                  values_to = "NES")
  
  # Order pathways by NES in first comparison
  pathway_order <- comparison_df$pathway_clean[order(comparison_df$NES1)]
  comparison_long$pathway_clean <- factor(comparison_long$pathway_clean, 
                                          levels = pathway_order)
  
  # Create grouped barplot
  p <- ggplot(comparison_long, aes(x = pathway_clean, y = NES, fill = Comparison)) +
    geom_col(position = "dodge", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("#C62828", "#1976D2"), 
                     name = "") +
    labs(title = paste0("Top ", top_n, " Pathways: ", name1, " vs ", name2),
         x = "",
         y = "Normalized Enrichment Score (NES)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "top",
          axis.text.y = element_text(size = 10))
  
  if(!is.null(filename)) {
    ggsave(filename, p, width = 12, height = 10)
  }
  
  print(p)
  return(p)
}

# JDM vs INFa comparison
cat("1. Creating JDM vs INFa comparison barplot...\n")
create_comparison_barplot(
  gsea_results$JDM_vs_Healthy, 
  gsea_results$INFa_vs_Control,
  "JDM (Disease)", "IFNa (Treatment)",
  top_n = 15,
  filename = "TenContrasts_GSEA/Barplot_JDM_vs_INFa_Top_Pathways.pdf"
)

# JDM vs INFb comparison
cat("2. Creating JDM vs INFb comparison barplot...\n")
create_comparison_barplot(
  gsea_results$JDM_vs_Healthy, 
  gsea_results$INFb_vs_Control,
  "JDM (Disease)", "IFNb (Treatment)",
  top_n = 15,
  filename = "TenContrasts_GSEA/Barplot_JDM_vs_INFb_Top_Pathways.pdf"
)

# INFa vs INFb comparison
cat("3. Creating INFa vs INFb comparison barplot...\n")
create_comparison_barplot(
  gsea_results$INFa_vs_Control, 
  gsea_results$INFb_vs_Control,
  "IFNa (Treatment)", "IFNb (Treatment)",
  top_n = 15,
  filename = "TenContrasts_GSEA/Barplot_INFa_vs_INFb_Top_Pathways.pdf"
)

cat("Comparison barplots saved\n\n")

################################################################################
# STEP 10: THREE-WAY COMPARISON BARPLOT
################################################################################

cat("Creating three-way comparison barplot...\n")

# Get all significant pathways
all_sig <- unique(c(jdm_pathways, ifna_pathways, ifnb_pathways))

if(length(all_sig) >= 10) {
  # Sort by absolute NES in JDM and take top 20
  jdm_order <- gsea_results$JDM_vs_Healthy[order(-abs(gsea_results$JDM_vs_Healthy$NES)), ]
  top_20_pathways <- head(jdm_order$pathway[jdm_order$pathway %in% all_sig], 20)
  
  # Create comparison data frame
  three_way_df <- data.frame(
    pathway = top_20_pathways,
    NES_JDM = gsea_results$JDM_vs_Healthy$NES[match(top_20_pathways, 
                                                     gsea_results$JDM_vs_Healthy$pathway)],
    NES_INFa = gsea_results$INFa_vs_Control$NES[match(top_20_pathways, 
                                                       gsea_results$INFa_vs_Control$pathway)],
    NES_INFb = gsea_results$INFb_vs_Control$NES[match(top_20_pathways, 
                                                       gsea_results$INFb_vs_Control$pathway)]
  )
  
  # Clean pathway names
  three_way_df$pathway_clean <- gsub("HALLMARK_", "", three_way_df$pathway)
  three_way_df$pathway_clean <- gsub("_", " ", three_way_df$pathway_clean)
  
  # Reshape for plotting
  three_way_long <- pivot_longer(three_way_df,
                                 cols = c(NES_JDM, NES_INFa, NES_INFb),
                                 names_to = "Comparison",
                                 values_to = "NES")
  
  three_way_long$Comparison <- gsub("NES_", "", three_way_long$Comparison)
  three_way_long$Comparison <- factor(three_way_long$Comparison, 
                                      levels = c("JDM", "INFa", "INFb"))
  
  # Order pathways by JDM NES
  pathway_order <- three_way_df$pathway_clean[order(three_way_df$NES_JDM)]
  three_way_long$pathway_clean <- factor(three_way_long$pathway_clean, 
                                          levels = pathway_order)
  
  # Create three-way barplot
  p_three_way <- ggplot(three_way_long, aes(x = pathway_clean, y = NES, fill = Comparison)) +
    geom_col(position = "dodge", width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c("JDM" = "#C62828", 
                                 "INFa" = "#1976D2", 
                                 "INFb" = "#F57C00"),
                     name = "Comparison",
                     labels = c("JDM (Disease)", 
                               "IFNa (Treatment)", 
                               "IFNb (Treatment)")) +
    labs(title = "Top 20 Pathways: Three-Way Comparison",
         subtitle = "JDM Disease Signature vs IFN Treatments",
         x = "",
         y = "Normalized Enrichment Score (NES)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.position = "top",
          axis.text.y = element_text(size = 9))
  
  print(p_three_way)
  ggsave("TenContrasts_GSEA/Barplot_ThreeWay_JDM_INFa_INFb_Top_Pathways.pdf", 
         p_three_way, width = 14, height = 12)
  
  cat("Three-way comparison barplot saved\n\n")
}

################################################################################
# STEP 11: HEATMAP OF NES SCORES
################################################################################

cat("Creating pathway enrichment heatmap...\n")

# Focus on key comparisons for heatmap
key_comparisons <- c("JDM_vs_Healthy", "INFa_vs_Control", "INFb_vs_Control",
                     "INFa_vs_JDM", "INFb_vs_JDM")

# Get all significant pathways from key comparisons
all_sig_key <- unique(unlist(sig_pathways_list[key_comparisons]))

if(length(all_sig_key) >= 2) {
  # Create matrix
  heatmap_data <- matrix(NA, nrow = length(all_sig_key), ncol = length(key_comparisons))
  rownames(heatmap_data) <- all_sig_key
  colnames(heatmap_data) <- key_comparisons
  
  # Fill in NES values
  for(i in seq_along(key_comparisons)) {
    comp <- key_comparisons[i]
    heatmap_data[, i] <- gsea_results[[comp]]$NES[match(all_sig_key, 
                                                         gsea_results[[comp]]$pathway)]
  }
  
  # Clean row names
  rownames(heatmap_data) <- gsub("HALLMARK_", "", rownames(heatmap_data))
  rownames(heatmap_data) <- gsub("_", " ", rownames(heatmap_data))
  
  # Clean column names
  colnames(heatmap_data) <- gsub("_", "\n", colnames(heatmap_data))
  
  # Create heatmap
  pdf("TenContrasts_GSEA/Heatmap_Key_Comparisons_Pathways.pdf", width = 10, height = 14)
  pheatmap(heatmap_data,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-3, 3, length.out = 101),
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           main = "Pathway Enrichment: Key Comparisons",
           fontsize_row = 8,
           fontsize_col = 10)
  dev.off()
  
  cat("Pathway heatmap saved\n\n")
}

################################################################################
# STEP 12: CREATE SUMMARY TABLE
################################################################################

pathway_summary <- data.frame(
  Comparison = names(gsea_results),
  Total_Pathways = sapply(gsea_results, nrow),
  Sig_Pathways = sapply(sig_pathways_list, length),
  Up_Pathways = sapply(sig_up_pathways_list, length),
  Down_Pathways = sapply(sig_down_pathways_list, length)
)

# Add recapitulation metrics
pathway_summary$JDM_Overlap <- NA
pathway_summary$JDM_Overlap_Pct <- NA

for(i in seq_along(pathway_summary$Comparison)) {
  comp <- as.character(pathway_summary$Comparison[i])
  if(comp != "JDM_vs_Healthy") {
    overlap <- length(intersect(jdm_pathways, sig_pathways_list[[comp]]))
    pathway_summary$JDM_Overlap[i] <- overlap
    pathway_summary$JDM_Overlap_Pct[i] <- round(100 * overlap / max(length(jdm_pathways), 1), 1)
  }
}

print(pathway_summary)
write.csv(pathway_summary, 
          "TenContrasts_GSEA/GSEA_Summary_All_Comparisons.csv", 
          row.names = FALSE)

# Save consensus pathways
if(length(overlap_jdm_ifna) > 0) {
  write.table(overlap_jdm_ifna, 
              "TenContrasts_GSEA/Consensus_Pathways_JDM_INFa.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

if(length(overlap_jdm_ifnb) > 0) {
  write.table(overlap_jdm_ifnb, 
              "TenContrasts_GSEA/Consensus_Pathways_JDM_INFb.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

################################################################################
# COMPLETE
################################################################################

cat("\n============================================================\n")
cat("GSEA ANALYSIS COMPLETE!\n")
cat("============================================================\n\n")

cat("Generated files in TenContrasts_GSEA/:\n")
cat("- 10 GSEA results CSV files\n")
cat("- 2 Venn diagrams (all pathways, upregulated)\n")
cat("- 2 NES correlation scatter plots\n")
cat("- 3 side-by-side comparison barplots\n")
cat("- 1 three-way comparison barplot\n")
cat("- 1 pathway enrichment heatmap\n")
cat("- 1 summary table\n")
cat("- Up to 2 consensus pathway lists\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("- JDM disease signature: %d significant pathways\n", length(jdm_pathways)))
cat(sprintf("- INFa treatment: %d significant pathways\n", length(ifna_pathways)))
cat(sprintf("- INFb treatment: %d significant pathways\n", length(ifnb_pathways)))
cat(sprintf("- JDM-INFa pathway overlap: %d (%.1f%% of JDM)\n", 
            length(overlap_jdm_ifna),
            100 * length(overlap_jdm_ifna) / max(length(jdm_pathways), 1)))
cat(sprintf("- JDM-INFb pathway overlap: %d (%.1f%% of JDM)\n", 
            length(overlap_jdm_ifnb),
            100 * length(overlap_jdm_ifnb) / max(length(jdm_pathways), 1)))
cat(sprintf("- NES correlation (JDM vs INFa): r = %.3f\n", cor_jdm_ifna))
cat(sprintf("- NES correlation (JDM vs INFb): r = %.3f\n\n", cor_jdm_ifnb))

cat("INTERPRETATION:\n")
if(length(overlap_jdm_ifnb) > length(overlap_jdm_ifna)) {
  cat("✓ IFNb shows better pathway-level recapitulation of JDM disease signature\n")
} else if(length(overlap_jdm_ifna) > length(overlap_jdm_ifnb)) {
  cat("✓ IFNa shows better pathway-level recapitulation of JDM disease signature\n")
} else {
  cat("✓ IFNa and IFNb show similar pathway-level recapitulation\n")
}

if(cor_jdm_ifnb > cor_jdm_ifna) {
  cat("✓ IFNb shows stronger correlation with JDM pathway enrichment patterns\n")
} else {
  cat("✓ IFNa shows stronger correlation with JDM pathway enrichment patterns\n")
}
