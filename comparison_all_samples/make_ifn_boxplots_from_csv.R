################################################################################
# Classic IFN Gene Boxplots - From IFN_genes_all_contrasts_LONG.csv
# X-axis: All contrasts (in custom order)
# Y-axis: log2FoldChange (or padj)
# Points: Individual IFN genes
#
# IMPORTANT: Control_vs_JDM is sign-flipped to represent JDM vs NT
################################################################################

library(ggplot2)
library(dplyr)

################################################################################
# CONFIGURATION
################################################################################

# Input file path
input_file <- "~/Documents/Covert_11358/all/IFN_Genes_Results/IFN_genes_all_contrasts_LONG.csv"

# Output directory
output_dir <- "IFN_Boxplots"
dir.create(output_dir, showWarnings = FALSE)

################################################################################
# 1. LOAD DATA
################################################################################

cat("Loading IFN gene data from:\n", input_file, "\n\n")

ifn_data <- read.csv(input_file, stringsAsFactors = FALSE)

cat("Data loaded successfully!\n")
cat("Total rows:", nrow(ifn_data), "\n")
cat("Unique genes:", length(unique(ifn_data$gene_symbol)), "\n")
cat("Unique contrasts:", length(unique(ifn_data$contrast)), "\n\n")

cat("Contrasts found:\n")
print(unique(ifn_data$contrast))
cat("\n")

################################################################################
# 2. PREPARE DATA
################################################################################

cat("Preparing data...\n")

# Remove rows with NA in key columns
ifn_plot_data <- ifn_data %>%
  filter(!is.na(log2FoldChange), !is.na(padj))

cat("Rows after removing NAs:", nrow(ifn_plot_data), "\n")

# Check distribution across contrasts
cat("\nData points per contrast:\n")
contrast_counts <- table(ifn_plot_data$contrast)
print(contrast_counts)
cat("\n")

################################################################################
# FLIP SIGN FOR Control_vs_JDM (to make it JDM vs NT)
################################################################################

cat("Flipping sign for Control_vs_JDM to make it JDM vs NT...\n")

ifn_plot_data <- ifn_plot_data %>%
  mutate(
    log2FoldChange = ifelse(contrast == "Control_vs_JDM", 
                            -log2FoldChange, 
                            log2FoldChange)
  )

cat("Sign flipped for Control_vs_JDM\n\n")

################################################################################
# CREATE CUSTOM LABELS AND ORDER
################################################################################

# Define mapping from original names to clean labels
contrast_labels <- c(
  "JDM_vs_Healthy" = "JDM vs Healthy",
  "Control_vs_Healthy" = "NT vs Healthy",
  "Control_vs_JDM" = "JDM vs NT",
  "INFa_vs_Control" = "INFa vs NT",
  "INFa_vs_Healthy" = "INFa vs Healthy",
  "INFa_vs_JDM" = "INFa vs JDM",
  "INFb_vs_Control" = "INFb vs NT",
  "INFb_vs_Healthy" = "INFb vs Healthy",
  "INFb_vs_JDM" = "INFb vs JDM",
  "INFb_vs_INFa" = "INFb vs INFa"
)

# Apply labels
ifn_plot_data$contrast_label <- contrast_labels[ifn_plot_data$contrast]

# Define custom order (as specified)
custom_order <- c(
  "JDM vs Healthy",
  "NT vs Healthy",
  "JDM vs NT",
  "INFa vs NT",
  "INFa vs Healthy",
  "INFa vs JDM",
  "INFb vs NT",
  "INFb vs Healthy",
  "INFb vs JDM",
  "INFb vs INFa"
)

# Set factor levels in custom order
ifn_plot_data$contrast_label <- factor(ifn_plot_data$contrast_label, 
                                       levels = custom_order)

cat("Applied custom labels and ordering\n\n")

# Calculate -log10(padj) for significance plot
ifn_plot_data$neg_log10_padj <- -log10(ifn_plot_data$padj)

# Cap infinite values (from padj = 0)
ifn_plot_data$neg_log10_padj[is.infinite(ifn_plot_data$neg_log10_padj)] <- 50

################################################################################
# 3. BOXPLOT 1: LOG2FOLDCHANGE BY CONTRAST
################################################################################

cat("Creating log2FoldChange boxplot...\n")

# Calculate plot width (1.2 inches per contrast, minimum 10)
n_contrasts <- length(unique(ifn_plot_data$contrast_label))
plot_width <- max(10, n_contrasts * 1.2)

cat("Number of contrasts:", n_contrasts, "\n")
cat("Plot width:", plot_width, "inches\n\n")

p1 <- ggplot(ifn_plot_data, aes(x = contrast_label, y = log2FoldChange)) +
  # Boxplot - white fill, black outline
  geom_boxplot(
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.75,
    width = 0.75
  ) +
  # Individual gene points
  geom_jitter(
    width = 0.15,
    alpha = 0.6,
    size = 4,
    color = "black",
    shape = 16
  ) +
  # Zero reference line
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.75
  ) +
  # Labels
  labs(
    x = NULL,
    y = "log2(Fold Change)"
  ) +
  # Classic theme
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 20,
      color = "black"
    ),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5,
      color = "black"
    ),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Save PDF
ggsave(
  file.path(output_dir, "IFN_genes_log2FC_boxplot.pdf"),
  plot = p1,
  width = plot_width,
  height = 6,
  bg = "white"
)

# Save PNG
ggsave(
  file.path(output_dir, "IFN_genes_log2FC_boxplot.png"),
  plot = p1,
  width = plot_width,
  height = 6,
  dpi = 300,
  bg = "white"
)

print(p1)

cat("Saved log2FC boxplot\n\n")

################################################################################
# 4. BOXPLOT 2: -LOG10(PADJ) BY CONTRAST
################################################################################

cat("Creating -log10(padj) boxplot...\n")

p2 <- ggplot(ifn_plot_data, aes(x = contrast_label, y = neg_log10_padj)) +
  # Boxplot - white fill, black outline
  geom_boxplot(
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.75,
    width = 0.75
  ) +
  # Individual gene points
  geom_jitter(
    width = 0.15,
    alpha = 0.6,
    size = 4,
    color = "black",
    shape = 16
  ) +
  # Significance threshold line
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.5
  ) +
  # Labels
  labs(
    x = NULL,
    y = "-log10(adjusted p-value)"
  ) +
  # Classic theme
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 20,
      color = "black"
    ),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    plot.title = element_text(
      face = "bold",
      size = 20,
      hjust = 0.5,
      color = "black"
    ),
    axis.line = element_line(color = "black", linewidth = 0.75),
    axis.ticks = element_line(color = "black", linewidth = 0.75),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Save PDF
ggsave(
  file.path(output_dir, "IFN_genes_padj_boxplot.pdf"),
  plot = p2,
  width = plot_width,
  height = 6,
  bg = "white"
)

# Save PNG
ggsave(
  file.path(output_dir, "IFN_genes_padj_boxplot.png"),
  plot = p2,
  width = plot_width,
  height = 6,
  dpi = 300,
  bg = "white"
)

print(p2)

cat("Saved padj boxplot\n\n")

################################################################################
# 5. SUMMARY STATISTICS
################################################################################

cat("Calculating summary statistics by contrast...\n\n")

summary_stats <- ifn_plot_data %>%
  group_by(contrast, contrast_label) %>%
  summarise(
    n_genes = n(),
    median_log2FC = round(median(log2FoldChange, na.rm = TRUE), 3),
    mean_log2FC = round(mean(log2FoldChange, na.rm = TRUE), 3),
    sd_log2FC = round(sd(log2FoldChange, na.rm = TRUE), 3),
    n_significant = sum(padj < 0.05, na.rm = TRUE),
    pct_significant = round(100 * n_significant / n_genes, 1),
    median_padj = signif(median(padj, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  # Keep the order from the factor levels
  arrange(factor(contrast_label, levels = custom_order))

print(summary_stats)

# Save summary table
write.csv(
  summary_stats,
  file.path(output_dir, "IFN_genes_summary_by_contrast.csv"),
  row.names = FALSE
)

################################################################################
# DONE!
################################################################################

cat("\n============================================================\n")
cat("BOXPLOTS COMPLETE!\n")
cat("============================================================\n\n")

cat("Summary:\n")
cat("- IFN genes analyzed:", length(unique(ifn_plot_data$gene_symbol)), "\n")
cat("- Contrasts plotted:", n_contrasts, "\n")
cat("- Total data points:", nrow(ifn_plot_data), "\n\n")

cat("Output files in", output_dir, "/:\n")
cat("1. IFN_genes_log2FC_boxplot.pdf\n")
cat("2. IFN_genes_log2FC_boxplot.png (300 dpi)\n")
cat("3. IFN_genes_padj_boxplot.pdf\n")
cat("4. IFN_genes_padj_boxplot.png (300 dpi)\n")
cat("5. IFN_genes_summary_by_contrast.csv\n\n")

cat("Contrasts in order (with median IFN gene log2FC):\n")
print(summary_stats[, c("contrast_label", "median_log2FC", "n_significant", "pct_significant")])
