library(tidyverse)
library(ggplot2)

# Optional package for non-overlapping gene labels on volcano plot
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# -----------------------------
# Load processed data
# -----------------------------
expr <- read.csv("data/processed/expression_filtered.csv", check.names = FALSE)
metadata <- read.csv("data/processed/metadata_matched.csv", check.names = FALSE)

# Set gene column as row names
rownames(expr) <- expr$gene
expr <- expr[, -1]

# Transpose so rows = samples, columns = genes
expr_t <- t(expr)

# Convert to numeric dataframe
expr_t <- as.data.frame(expr_t)
expr_t[] <- lapply(expr_t, as.numeric)

# Remove zero-variance genes
gene_variance <- apply(expr_t, 2, var, na.rm = TRUE)
expr_t <- expr_t[, gene_variance > 0]

cat("Number of genes after removing zero-variance genes:", ncol(expr_t), "\n")

# -----------------------------
# PCA Analysis
# -----------------------------
pca <- prcomp(expr_t, scale. = TRUE)

# Variance explained
pca_var <- pca$sdev^2
pca_var_explained <- round(100 * pca_var / sum(pca_var), 1)

# Create PCA dataframe
pca_df <- as.data.frame(pca$x)
pca_df$sample_id <- rownames(pca_df)

# Standardise sample IDs
pca_df$sample_id <- gsub("\\.", "-", pca_df$sample_id)
pca_df$sample_id <- gsub("^X", "", pca_df$sample_id)
pca_df$sample_id <- substr(pca_df$sample_id, 1, 15)

metadata$sample_id <- gsub("\\.", "-", metadata$sample_id)
metadata$sample_id <- gsub("^X", "", metadata$sample_id)
metadata$sample_id <- substr(metadata$sample_id, 1, 15)

# Debug checks
cat("First 5 PCA sample IDs:\n")
print(head(pca_df$sample_id, 5))

cat("First 5 metadata sample IDs:\n")
print(head(metadata$sample_id, 5))

cat("Number of overlapping sample IDs:\n")
print(length(intersect(pca_df$sample_id, metadata$sample_id)))

# Merge PCA coordinates with metadata
pca_df <- pca_df %>%
  left_join(metadata, by = "sample_id")

# Check merge worked
cat("Sample type counts after merge:\n")
print(table(pca_df$sample_type, useNA = "ifany"))

# Improved PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = sample_type)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = c("Normal" = "#E64B35", "Tumor" = "#00A087")) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA of TCGA LUAD Expression Data",
    x = paste0("PC1 (", pca_var_explained[1], "%)"),
    y = paste0("PC2 (", pca_var_explained[2], "%)"),
    color = "Sample Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

print(pca_plot)

# Save PCA plot
ggsave("results/figures/pca_plot.png", plot = pca_plot, width = 7, height = 5)

# -----------------------------
# Final Polished Volcano Plot
# -----------------------------

de_results <- read.csv("results/tables/differential_expression_results.csv", check.names = FALSE)

de_results <- de_results %>%
  mutate(
    neg_log10_adjP = -log10(adj.P.Val),
    significance = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# Top genes (more controlled selection)
top_genes <- de_results %>%
  filter(adj.P.Val < 1e-20) %>%
  arrange(adj.P.Val) %>%
  head(8)

volcano_plot <- ggplot(de_results, aes(x = logFC, y = neg_log10_adjP, color = significance)) +
  geom_point(alpha = 0.7, size = 1.3) +

  # Threshold lines
geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.6) +
geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.6) +

  # Colours
  scale_color_manual(values = c(
    "Upregulated" = "#3C8DFF",
    "Downregulated" = "#E64B35",
    "Not Significant" = "grey75"
  )) +

  theme_classic(base_size = 15) +

  labs(
    title = "Differential Expression in TCGA LUAD",
    subtitle = "Tumour vs Normal",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P Value",
    color = "Gene Category"
  ) +

  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

# Add labels if ggrepel exists
if (has_ggrepel) {
  volcano_plot <- volcano_plot +
    ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.4,
      point.padding = 0.3
    )
}

print(volcano_plot)

ggsave(
  "results/figures/volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6
)

# -----------------------------
# Heatmap of Top Differentially Expressed Genes
# -----------------------------

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Package 'pheatmap' is not installed. Run install.packages('pheatmap')")
}
library(pheatmap)

cat("\nStarting heatmap generation...\n")

# Load full expression matrix again
expr_heatmap <- read.csv("data/processed/expression_filtered.csv", check.names = FALSE)
metadata_heatmap <- read.csv("data/processed/metadata_matched.csv", check.names = FALSE)

# Set gene names as row names
rownames(expr_heatmap) <- expr_heatmap$gene
expr_heatmap <- expr_heatmap[, -1]

# Standardise sample IDs
colnames(expr_heatmap) <- gsub("\\.", "-", colnames(expr_heatmap))
colnames(expr_heatmap) <- gsub("^X", "", colnames(expr_heatmap))
colnames(expr_heatmap) <- substr(colnames(expr_heatmap), 1, 15)

metadata_heatmap$sample_id <- gsub("\\.", "-", metadata_heatmap$sample_id)
metadata_heatmap$sample_id <- gsub("^X", "", metadata_heatmap$sample_id)
metadata_heatmap$sample_id <- substr(metadata_heatmap$sample_id, 1, 15)

# Keep shared samples and match order
common_samples_hm <- intersect(colnames(expr_heatmap), metadata_heatmap$sample_id)
cat("Heatmap shared samples:", length(common_samples_hm), "\n")

expr_heatmap <- expr_heatmap[, common_samples_hm, drop = FALSE]
metadata_heatmap <- metadata_heatmap %>%
  filter(sample_id %in% common_samples_hm) %>%
  slice(match(colnames(expr_heatmap), sample_id))

# Select top genes from differential expression results
top_heatmap_genes <- de_results %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 2) %>%   # stronger genes only
  arrange(adj.P.Val) %>%
  head(30) %>%                                    # fewer = cleaner heatmap
  pull(gene)

cat("Top heatmap genes before intersect:", length(top_heatmap_genes), "\n")

# Keep only genes present in expression matrix
top_heatmap_genes <- intersect(top_heatmap_genes, rownames(expr_heatmap))
cat("Top heatmap genes after intersect:", length(top_heatmap_genes), "\n")

if (length(top_heatmap_genes) == 0) {
  stop("No top heatmap genes found in expression matrix.")
}

# Subset expression matrix
heatmap_matrix <- expr_heatmap[top_heatmap_genes, , drop = FALSE]

# Convert to numeric matrix
heatmap_matrix <- as.matrix(heatmap_matrix)
mode(heatmap_matrix) <- "numeric"

cat("Heatmap matrix dimensions:", dim(heatmap_matrix), "\n")

# Remove genes with zero variance before scaling
row_vars <- apply(heatmap_matrix, 1, var, na.rm = TRUE)
heatmap_matrix <- heatmap_matrix[row_vars > 0, , drop = FALSE]

cat("Heatmap matrix dimensions after removing zero-variance rows:", dim(heatmap_matrix), "\n")

if (nrow(heatmap_matrix) == 0) {
  stop("No variable genes left for heatmap after filtering zero-variance rows.")
}

# Scale by gene (row)
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

# Replace any non-finite values
heatmap_matrix_scaled[!is.finite(heatmap_matrix_scaled)] <- 0

# Annotation for columns (samples)
annotation_col <- data.frame(
  SampleType = metadata_heatmap$sample_type
)
rownames(annotation_col) <- metadata_heatmap$sample_id

# Colour mapping for annotations
annotation_colors <- list(
  SampleType = c(
    "Normal" = "#D55E00",
    "Tumor" = "#009E73"
  )
)

# Output path
heatmap_file <- "results/figures/heatmap_top50_genes.png"
cat("Saving heatmap to:", heatmap_file, "\n")

pheatmap(
  heatmap_matrix_scaled,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,

  # cleaner display
  show_rownames = TRUE,
  show_colnames = FALSE,

  # clustering
  cluster_rows = TRUE,
  cluster_cols = TRUE,

  # better colours (publication style)
  color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),

  # aesthetics
  fontsize_row = 9,
  border_color = NA,

  # title
  main = "Top Differentially Expressed Genes (TCGA LUAD)",

  # save
  filename = heatmap_file,
  width = 10,
  height = 8
)

cat("Heatmap generation complete.\n")