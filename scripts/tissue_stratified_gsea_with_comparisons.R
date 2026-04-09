################################################################################
# GSEA Analysis for Tissue-Stratified Results
# Purpose: Run GSEA on JDMS, IFNa, and IFNb results and create comparison plots
# Style: Barplot comparisons
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

cat("=== TISSUE-STRATIFIED GSEA ANALYSIS ===\n\n")

cat("Loading Hallmark gene sets from MSigDB...\n")
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- hallmark_sets %>% 
  split(x = .$gene_symbol, f = .$gs_name)

cat("Loaded", length(hallmark_list), "Hallmark gene sets\n\n")

################################################################################
# STEP 2: LOAD TISSUE-STRATIFIED RESULTS (WITH GENE SYMBOLS)
################################################################################

cat("Loading tissue-stratified DESeq2 results...\n")

# These files should exist after running complete_tissue_stratified_analysis_with_symbols.R
res_jdms <- read.csv("DESeq2_PatientMuscle_JDMS_vs_Healthy.csv")
res_ifna <- read.csv("DESeq2_Myobundles_IFNa_vs_NT.csv")
res_ifnb <- read.csv("DESeq2_Myobundles_IFNb_vs_NT.csv")

cat("JDMS results:", nrow(res_jdms), "genes\n")
cat("IFNa results:", nrow(res_ifna), "genes\n")
cat("IFNb results:", nrow(res_ifnb), "genes\n\n")

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
  
  # Calculate ranking metric: sign(log2FC) * -log10(pvalue)
  results_df$rank_metric <- sign(results_df$log2FoldChange) * -log10(results_df$pvalue)
  
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

# Create ranked lists
ranked_jdms <- prepare_ranked_list(res_jdms)
ranked_ifna <- prepare_ranked_list(res_ifna)
ranked_ifnb <- prepare_ranked_list(res_ifnb)

cat("JDMS ranked genes:", length(ranked_jdms), "\n")
cat("IFNa ranked genes:", length(ranked_ifna), "\n")
cat("IFNb ranked genes:", length(ranked_ifnb), "\n\n")

################################################################################
# STEP 4: RUN FGSEA
################################################################################

cat("Running GSEA on all comparisons...\n\n")

# JDMS vs Healthy (Patient Muscle)
cat("1. Running GSEA: JDMS vs Healthy (Patient Muscle)...\n")
gsea_jdms <- fgsea(
  pathways = hallmark_list,
  stats = ranked_jdms,
  minSize = 15,
  maxSize = 500
)
gsea_jdms <- gsea_jdms[order(gsea_jdms$NES, decreasing = TRUE), ]

sig_up_jdms <- sum(gsea_jdms$padj < 0.05 & gsea_jdms$NES > 0, na.rm = TRUE)
sig_down_jdms <- sum(gsea_jdms$padj < 0.05 & gsea_jdms$NES < 0, na.rm = TRUE)
cat("   Significant pathways:", sig_up_jdms + sig_down_jdms, 
    "(", sig_up_jdms, "up,", sig_down_jdms, "down)\n\n")

# IFNa vs NT (Myobundles)
cat("2. Running GSEA: IFNa vs NT (Myobundles)...\n")
gsea_ifna <- fgsea(
  pathways = hallmark_list,
  stats = ranked_ifna,
  minSize = 15,
  maxSize = 500
)
gsea_ifna <- gsea_ifna[order(gsea_ifna$NES, decreasing = TRUE), ]

sig_up_ifna <- sum(gsea_ifna$padj < 0.05 & gsea_ifna$NES > 0, na.rm = TRUE)
sig_down_ifna <- sum(gsea_ifna$padj < 0.05 & gsea_ifna$NES < 0, na.rm = TRUE)
cat("   Significant pathways:", sig_up_ifna + sig_down_ifna, 
    "(", sig_up_ifna, "up,", sig_down_ifna, "down)\n\n")

# IFNb vs NT (Myobundles)
cat("3. Running GSEA: IFNb vs NT (Myobundles)...\n")
gsea_ifnb <- fgsea(
  pathways = hallmark_list,
  stats = ranked_ifnb,
  minSize = 15,
  maxSize = 500
)
gsea_ifnb <- gsea_ifnb[order(gsea_ifnb$NES, decreasing = TRUE), ]

sig_up_ifnb <- sum(gsea_ifnb$padj < 0.05 & gsea_ifnb$NES > 0, na.rm = TRUE)
sig_down_ifnb <- sum(gsea_ifnb$padj < 0.05 & gsea_ifnb$NES < 0, na.rm = TRUE)
cat("   Significant pathways:", sig_up_ifnb + sig_down_ifnb, 
    "(", sig_up_ifnb, "up,", sig_down_ifnb, "down)\n\n")

################################################################################
# STEP 5: SAVE GSEA RESULTS
################################################################################

cat("Saving GSEA results to CSV files...\n")

# Convert leadingEdge lists to strings
gsea_jdms_df <- as.data.frame(gsea_jdms)
gsea_jdms_df$leadingEdge <- sapply(gsea_jdms_df$leadingEdge, 
                                   function(x) paste(x, collapse = ";"))
write.csv(gsea_jdms_df, "GSEA_PatientMuscle_JDMS_vs_Healthy.csv", row.names = FALSE)

gsea_ifna_df <- as.data.frame(gsea_ifna)
gsea_ifna_df$leadingEdge <- sapply(gsea_ifna_df$leadingEdge, 
                                   function(x) paste(x, collapse = ";"))
write.csv(gsea_ifna_df, "GSEA_Myobundles_IFNa_vs_NT.csv", row.names = FALSE)

gsea_ifnb_df <- as.data.frame(gsea_ifnb)
gsea_ifnb_df$leadingEdge <- sapply(gsea_ifnb_df$leadingEdge, 
                                   function(x) paste(x, collapse = ";"))
write.csv(gsea_ifnb_df, "GSEA_Myobundles_IFNb_vs_NT.csv", row.names = FALSE)

cat("GSEA results saved\n\n")

################################################################################
# STEP 6: PATHWAY OVERLAP ANALYSIS
################################################################################

cat("=== PATHWAY OVERLAP ANALYSIS ===\n\n")

# Get significant pathways from each comparison
sig_pathways_jdms <- gsea_jdms$pathway[gsea_jdms$padj < 0.05]
sig_pathways_ifna <- gsea_ifna$pathway[gsea_ifna$padj < 0.05]
sig_pathways_ifnb <- gsea_ifnb$pathway[gsea_ifnb$padj < 0.05]

# Get upregulated pathways
sig_up_pathways_jdms <- gsea_jdms$pathway[gsea_jdms$padj < 0.05 & gsea_jdms$NES > 0]
sig_up_pathways_ifna <- gsea_ifna$pathway[gsea_ifna$padj < 0.05 & gsea_ifna$NES > 0]
sig_up_pathways_ifnb <- gsea_ifnb$pathway[gsea_ifnb$padj < 0.05 & gsea_ifnb$NES > 0]

# Get downregulated pathways
sig_down_pathways_jdms <- gsea_jdms$pathway[gsea_jdms$padj < 0.05 & gsea_jdms$NES < 0]
sig_down_pathways_ifna <- gsea_ifna$pathway[gsea_ifna$padj < 0.05 & gsea_ifna$NES < 0]
sig_down_pathways_ifnb <- gsea_ifnb$pathway[gsea_ifnb$padj < 0.05 & gsea_ifnb$NES < 0]

# Calculate overlaps
overlap_jdms_ifna <- intersect(sig_pathways_jdms, sig_pathways_ifna)
overlap_jdms_ifnb <- intersect(sig_pathways_jdms, sig_pathways_ifnb)
overlap_ifna_ifnb <- intersect(sig_pathways_ifna, sig_pathways_ifnb)
overlap_all_three <- Reduce(intersect, list(sig_pathways_jdms, sig_pathways_ifna, sig_pathways_ifnb))

# Print overlap summary
cat("PATHWAY OVERLAP SUMMARY:\n")
cat("JDMS significant pathways:", length(sig_pathways_jdms), "\n")
cat("IFNa significant pathways:", length(sig_pathways_ifna), "\n")
cat("IFNb significant pathways:", length(sig_pathways_ifnb), "\n\n")

cat("JDMS-IFNa overlap:", length(overlap_jdms_ifna), 
    "(", round(100 * length(overlap_jdms_ifna) / length(sig_pathways_jdms), 1), "% of JDMS)\n")
cat("JDMS-IFNb overlap:", length(overlap_jdms_ifnb), 
    "(", round(100 * length(overlap_jdms_ifnb) / length(sig_pathways_jdms), 1), "% of JDMS)\n")
cat("IFNa-IFNb overlap:", length(overlap_ifna_ifnb), "\n")
cat("All three overlap:", length(overlap_all_three), "\n\n")

if(length(overlap_jdms_ifna) > 0) {
  cat("CONSENSUS PATHWAYS (JDMS-IFNa):\n")
  for(p in overlap_jdms_ifna) {
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

# 3-way Venn - All significant pathways
if(length(sig_pathways_jdms) > 0 && length(sig_pathways_ifna) > 0 && length(sig_pathways_ifnb) > 0) {
  venn.plot.all <- venn.diagram(
    x = list(JDMS = sig_pathways_jdms, 
             IFNa = sig_pathways_ifna, 
             IFNb = sig_pathways_ifnb),
    category.names = c("JDMS\n(Muscle)", "IFNa\n(Myobundle)", "IFNb\n(Myobundle)"),
    filename = NULL,
    output = TRUE,
    main = "Pathway Overlap: JDMS vs IFN Treatments",
    main.cex = 1.5,
    lwd = 2,
    col = c("#C62828", "#1976D2", "#F57C00"),
    fill = c(alpha("#C62828", 0.3), alpha("#1976D2", 0.3), alpha("#F57C00", 0.3)),
    cex = 1.5,
    fontface = "bold",
    cat.cex = 1.3,
    cat.fontface = "bold"
  )
  
  pdf("Venn_JDMS_vs_IFNa_vs_IFNb_Pathways_ALL.pdf", width = 10, height = 10)
  grid.draw(venn.plot.all)
  dev.off()
}

# 3-way Venn - Upregulated pathways
if(length(sig_up_pathways_jdms) > 0 && length(sig_up_pathways_ifna) > 0 && length(sig_up_pathways_ifnb) > 0) {
  venn.plot.up <- venn.diagram(
    x = list(JDMS = sig_up_pathways_jdms, 
             IFNa = sig_up_pathways_ifna, 
             IFNb = sig_up_pathways_ifnb),
    category.names = c("JDMS\n(Muscle)", "IFNa\n(Myobundle)", "IFNb\n(Myobundle)"),
    filename = NULL,
    output = TRUE,
    main = "Upregulated Pathways: JDMS vs IFN Treatments",
    main.cex = 1.5,
    lwd = 2,
    col = c("#C62828", "#1976D2", "#F57C00"),
    fill = c(alpha("#C62828", 0.3), alpha("#1976D2", 0.3), alpha("#F57C00", 0.3)),
    cex = 1.5,
    fontface = "bold",
    cat.cex = 1.3,
    cat.fontface = "bold"
  )
  
  pdf("Venn_JDMS_vs_IFNa_vs_IFNb_Pathways_UP.pdf", width = 10, height = 10)
  grid.draw(venn.plot.up)
  dev.off()
}

cat("Venn diagrams saved\n\n")

################################################################################
# STEP 8: CORRELATION OF NES SCORES
################################################################################

cat("Analyzing correlation of pathway enrichment scores...\n")

# Merge all three GSEA results
merged_gsea <- merge(
  gsea_jdms[, c("pathway", "NES", "padj")],
  gsea_ifna[, c("pathway", "NES", "padj")],
  by = "pathway",
  suffixes = c("_JDMS", "_IFNa")
)

merged_gsea <- merge(
  merged_gsea,
  gsea_ifnb[, c("pathway", "NES", "padj")],
  by = "pathway"
)
names(merged_gsea)[names(merged_gsea) == "NES"] <- "NES_IFNb"
names(merged_gsea)[names(merged_gsea) == "padj"] <- "padj_IFNb"

# Calculate correlations
cor_jdms_ifna <- cor(merged_gsea$NES_JDMS, merged_gsea$NES_IFNa)
cor_jdms_ifnb <- cor(merged_gsea$NES_JDMS, merged_gsea$NES_IFNb)
cor_ifna_ifnb <- cor(merged_gsea$NES_IFNa, merged_gsea$NES_IFNb)

cat("Correlation of NES (all pathways):\n")
cat("  JDMS vs IFNa:", round(cor_jdms_ifna, 3), "\n")
cat("  JDMS vs IFNb:", round(cor_jdms_ifnb, 3), "\n")
cat("  IFNa vs IFNb:", round(cor_ifna_ifnb, 3), "\n\n")

# Scatter plots
p_jdms_ifna <- ggplot(merged_gsea, aes(x = NES_JDMS, y = NES_IFNa)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Pathway Enrichment: JDMS vs IFNa",
       subtitle = paste0("Pearson r = ", round(cor_jdms_ifna, 3), ", p = ", format.pval(cor.test(merged_gsea$NES_JDMS, merged_gsea$NES_IFNa)$p.value, digits = 2)),
       x = "NES (JDMS - Patient Muscle)",
       y = "NES (IFNa - Myobundles)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("Scatter_JDMS_vs_IFNa_Pathway_NES.pdf", p_jdms_ifna, width = 8, height = 8)

p_jdms_ifnb <- ggplot(merged_gsea, aes(x = NES_JDMS, y = NES_IFNb)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Pathway Enrichment: JDMS vs IFNb",
       subtitle = paste0("Pearson r = ", round(cor_jdms_ifnb, 3), ", p = ", format.pval(cor.test(merged_gsea$NES_JDMS, merged_gsea$NES_IFNb)$p.value, digits = 2)),
       x = "NES (JDMS - Patient Muscle)",
       y = "NES (IFNb - Myobundles)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("Scatter_JDMS_vs_IFNb_Pathway_NES.pdf", p_jdms_ifnb, width = 8, height = 8)

cat("Correlation plots saved\n\n")

################################################################################
# STEP 9: SIDE-BY-SIDE COMPARISON BARPLOTS (DUKE VS UCSF STYLE)
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

# JDMS vs IFNa comparison (top JDMS pathways)
cat("1. Creating JDMS vs IFNa comparison barplot...\n")
create_comparison_barplot(
  gsea_jdms, gsea_ifna,
  "JDMS (Muscle)", "IFNa (Myobundle)",
  top_n = 15,
  filename = "Barplot_JDMS_vs_IFNa_Top_Pathways.pdf"
)

# JDMS vs IFNb comparison (top JDMS pathways)
cat("2. Creating JDMS vs IFNb comparison barplot...\n")
create_comparison_barplot(
  gsea_jdms, gsea_ifnb,
  "JDMS (Muscle)", "IFNb (Myobundle)",
  top_n = 15,
  filename = "Barplot_JDMS_vs_IFNb_Top_Pathways.pdf"
)

# IFNa vs IFNb comparison (top IFNa pathways)
cat("3. Creating IFNa vs IFNb comparison barplot...\n")
create_comparison_barplot(
  gsea_ifna, gsea_ifnb,
  "IFNa (Myobundle)", "IFNb (Myobundle)",
  top_n = 15,
  filename = "Barplot_IFNa_vs_IFNb_Top_Pathways.pdf"
)

cat("Comparison barplots saved\n\n")

################################################################################
# STEP 10: THREE-WAY COMPARISON BARPLOT (ALL COMPARISONS)
################################################################################

cat("Creating three-way comparison barplot...\n")

# Get pathways that are significant in at least one comparison
all_sig_pathways <- unique(c(sig_pathways_jdms, sig_pathways_ifna, sig_pathways_ifnb))

if(length(all_sig_pathways) >= 10) {
  # Sort by absolute NES in JDMS and take top 20
  jdms_order <- gsea_jdms[order(-abs(gsea_jdms$NES)), ]
  top_20_pathways <- head(jdms_order$pathway[jdms_order$pathway %in% all_sig_pathways], 20)
  
  # Create comparison data frame
  three_way_df <- data.frame(
    pathway = top_20_pathways,
    NES_JDMS = gsea_jdms$NES[match(top_20_pathways, gsea_jdms$pathway)],
    NES_IFNa = gsea_ifna$NES[match(top_20_pathways, gsea_ifna$pathway)],
    NES_IFNb = gsea_ifnb$NES[match(top_20_pathways, gsea_ifnb$pathway)]
  )
  
  # Clean pathway names
  three_way_df$pathway_clean <- gsub("HALLMARK_", "", three_way_df$pathway)
  three_way_df$pathway_clean <- gsub("_", " ", three_way_df$pathway_clean)
  
  # Reshape for plotting
  three_way_long <- pivot_longer(three_way_df,
                                 cols = c(NES_JDMS, NES_IFNa, NES_IFNb),
                                 names_to = "Comparison",
                                 values_to = "NES")
  
  three_way_long$Comparison <- gsub("NES_", "", three_way_long$Comparison)
  three_way_long$Comparison <- factor(three_way_long$Comparison, 
                                      levels = c("JDMS", "IFNa", "IFNb"))
  
  # Order pathways by JDMS NES
  pathway_order <- three_way_df$pathway_clean[order(three_way_df$NES_JDMS)]
  three_way_long$pathway_clean <- factor(three_way_long$pathway_clean, 
                                          levels = pathway_order)
  
  # Create three-way barplot
  p_three_way <- ggplot(three_way_long, aes(x = pathway_clean, y = NES, fill = Comparison)) +
    geom_col(position = "dodge", width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c("JDMS" = "#C62828", 
                                 "IFNa" = "#1976D2", 
                                 "IFNb" = "#F57C00"),
                     name = "Comparison",
                     labels = c("JDMS (Muscle)", 
                               "IFNa (Myobundle)", 
                               "IFNb (Myobundle)")) +
    labs(title = "Top 20 Pathways: Three-Way Comparison",
         subtitle = "JDMS (Patient Muscle) vs IFNa/IFNb (Myobundles)",
         x = "",
         y = "Normalized Enrichment Score (NES)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.position = "top",
          axis.text.y = element_text(size = 9))
  
  print(p_three_way)
  ggsave("Barplot_ThreeWay_JDMS_IFNa_IFNb_Top_Pathways.pdf", 
         p_three_way, width = 14, height = 12)
  
  cat("Three-way comparison barplot saved\n\n")
}

################################################################################
# STEP 11: HEATMAP OF NES SCORES
################################################################################

cat("Creating pathway enrichment heatmap...\n")

# Get all significant pathways from any comparison
all_sig_pathways <- unique(c(sig_pathways_jdms, sig_pathways_ifna, sig_pathways_ifnb))

if(length(all_sig_pathways) >= 2) {
  # Create matrix
  heatmap_data <- matrix(NA, nrow = length(all_sig_pathways), ncol = 3)
  rownames(heatmap_data) <- all_sig_pathways
  colnames(heatmap_data) <- c("JDMS\n(Muscle)", "IFNa\n(Myobundle)", "IFNb\n(Myobundle)")
  
  # Fill in NES values
  heatmap_data[, 1] <- gsea_jdms$NES[match(all_sig_pathways, gsea_jdms$pathway)]
  heatmap_data[, 2] <- gsea_ifna$NES[match(all_sig_pathways, gsea_ifna$pathway)]
  heatmap_data[, 3] <- gsea_ifnb$NES[match(all_sig_pathways, gsea_ifnb$pathway)]
  
  # Clean row names
  rownames(heatmap_data) <- gsub("HALLMARK_", "", rownames(heatmap_data))
  rownames(heatmap_data) <- gsub("_", " ", rownames(heatmap_data))
  
  # Create heatmap
  pdf("Heatmap_JDMS_IFNa_IFNb_Pathways.pdf", width = 8, height = 12)
  pheatmap(heatmap_data,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-3, 3, length.out = 101),
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           main = "Pathway Enrichment: JDMS vs IFN Treatments",
           fontsize_row = 8,
           fontsize_col = 10)
  dev.off()
  
  cat("Pathway heatmap saved\n\n")
}

################################################################################
# STEP 12: CREATE SUMMARY TABLE
################################################################################

pathway_summary <- data.frame(
  Metric = c("Significant pathways (JDMS)",
             "Significant pathways (IFNa)",
             "Significant pathways (IFNb)",
             "JDMS-IFNa overlap",
             "JDMS-IFNb overlap",
             "IFNa-IFNb overlap",
             "All three overlap",
             "JDMS-IFNa overlap % (of JDMS)",
             "JDMS-IFNb overlap % (of JDMS)",
             "Upregulated pathways (JDMS)",
             "Upregulated pathways (IFNa)",
             "Upregulated pathways (IFNb)",
             "Downregulated pathways (JDMS)",
             "Downregulated pathways (IFNa)",
             "Downregulated pathways (IFNb)",
             "NES correlation (JDMS vs IFNa)",
             "NES correlation (JDMS vs IFNb)",
             "NES correlation (IFNa vs IFNb)"),
  Value = c(length(sig_pathways_jdms),
            length(sig_pathways_ifna),
            length(sig_pathways_ifnb),
            length(overlap_jdms_ifna),
            length(overlap_jdms_ifnb),
            length(overlap_ifna_ifnb),
            length(overlap_all_three),
            round(100 * length(overlap_jdms_ifna) / max(length(sig_pathways_jdms), 1), 1),
            round(100 * length(overlap_jdms_ifnb) / max(length(sig_pathways_jdms), 1), 1),
            length(sig_up_pathways_jdms),
            length(sig_up_pathways_ifna),
            length(sig_up_pathways_ifnb),
            length(sig_down_pathways_jdms),
            length(sig_down_pathways_ifna),
            length(sig_down_pathways_ifnb),
            round(cor_jdms_ifna, 3),
            round(cor_jdms_ifnb, 3),
            round(cor_ifna_ifnb, 3))
)

print(pathway_summary)
write.csv(pathway_summary, "GSEA_Pathway_Summary_JDMS_vs_IFN.csv", row.names = FALSE)

# Save consensus pathways
if(length(overlap_jdms_ifna) > 0) {
  write.table(overlap_jdms_ifna, "Consensus_Pathways_JDMS_IFNa.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

if(length(overlap_jdms_ifnb) > 0) {
  write.table(overlap_jdms_ifnb, "Consensus_Pathways_JDMS_IFNb.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

if(length(overlap_all_three) > 0) {
  write.table(overlap_all_three, "Consensus_Pathways_All_Three.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

cat("\n========================================\n")
cat("TISSUE-STRATIFIED GSEA ANALYSIS COMPLETE!\n")
cat("========================================\n\n")

cat("Generated files:\n")
cat("- 3 GSEA results CSV files\n")
cat("- 2 Venn diagrams (all pathways, upregulated)\n")
cat("- 2 NES correlation scatter plots\n")
cat("- 3 side-by-side comparison barplots\n")
cat("- 1 three-way comparison barplot\n")
cat("- 1 pathway enrichment heatmap\n")
cat("- 1 pathway summary table\n")
cat("- Up to 3 consensus pathway lists\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("- %d%% of JDMS pathways overlap with IFNa\n", 
            round(100 * length(overlap_jdms_ifna) / max(length(sig_pathways_jdms), 1), 1)))
cat(sprintf("- %d%% of JDMS pathways overlap with IFNb\n", 
            round(100 * length(overlap_jdms_ifnb) / max(length(sig_pathways_jdms), 1), 1)))
cat(sprintf("- NES correlation (JDMS vs IFNa): r = %.3f\n", cor_jdms_ifna))
cat(sprintf("- NES correlation (JDMS vs IFNb): r = %.3f\n", cor_jdms_ifnb))
