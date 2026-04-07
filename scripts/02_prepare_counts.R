library(tidyverse)

# Load expression data
expr <- read.delim("data/raw/expression/tcga_luad_expression.tsv", check.names = FALSE)

# Load metadata
metadata <- read.csv("data/processed/metadata_clean.csv")

# Rename first column
colnames(expr)[1] <- "gene"

# Keep only Tumor and Normal
metadata_filtered <- metadata %>%
  filter(sample_type %in% c("Tumor", "Normal"))

# Extract sample IDs
expr_samples <- colnames(expr)[-1]

# Check overlap
common_samples <- intersect(expr_samples, metadata_filtered$sample_id)

cat("Total expression samples:", length(expr_samples), "\n")
cat("Metadata samples:", nrow(metadata_filtered), "\n")
cat("Common samples:", length(common_samples), "\n")

# Filter metadata
metadata_matched <- metadata_filtered %>%
  filter(sample_id %in% common_samples)

# Reorder metadata to match expression
metadata_matched <- metadata_matched %>%
  arrange(match(sample_id, expr_samples))

# Filter expression matrix
expr_filtered <- expr %>%
  select(gene, all_of(metadata_matched$sample_id))

# Final sanity check
stopifnot(all(colnames(expr_filtered)[-1] == metadata_matched$sample_id))

# Save outputs
write.csv(expr_filtered, "data/processed/expression_filtered.csv", row.names = FALSE)
write.csv(metadata_matched, "data/processed/metadata_matched.csv", row.names = FALSE)

# Print summary
cat("Final expression dimensions:", dim(expr_filtered), "\n")
cat("Final metadata dimensions:", dim(metadata_matched), "\n")

# Preview
head(expr_filtered[, 1:6])
head(metadata_matched)