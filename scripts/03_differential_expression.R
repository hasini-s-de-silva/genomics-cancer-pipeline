library(tidyverse)
library(limma)

# Load processed data
expr <- read.csv("data/processed/expression_filtered.csv", check.names = FALSE)
metadata <- read.csv("data/processed/metadata_matched.csv", check.names = FALSE)

# Set gene column as row names
rownames(expr) <- expr$gene
expr <- expr[, -1]

# Standardise sample IDs
metadata$sample_id <- gsub("\\.", "-", metadata$sample_id)
metadata$sample_id <- gsub("^X", "", metadata$sample_id)
metadata$sample_id <- substr(metadata$sample_id, 1, 15)

colnames(expr) <- gsub("\\.", "-", colnames(expr))
colnames(expr) <- gsub("^X", "", colnames(expr))
colnames(expr) <- substr(colnames(expr), 1, 15)

# Keep only shared samples
common_samples <- intersect(colnames(expr), metadata$sample_id)

expr <- expr[, common_samples]
metadata <- metadata %>%
  filter(sample_id %in% common_samples)

# Reorder metadata to match expression columns
metadata <- metadata %>%
  slice(match(colnames(expr), sample_id))

# Final check
stopifnot(all(colnames(expr) == metadata$sample_id))

# Create group factor
metadata$sample_type <- factor(metadata$sample_type, levels = c("Normal", "Tumor"))

# Build design matrix
design <- model.matrix(~ sample_type, data = metadata)
colnames(design) <- c("Intercept", "Tumor_vs_Normal")

# Fit linear model
fit <- lmFit(as.matrix(expr), design)
fit <- eBayes(fit)

# Extract results
results <- topTable(
  fit,
  coef = "Tumor_vs_Normal",
  number = Inf,
  sort.by = "P"
)

# Add gene column
results <- results %>%
  rownames_to_column("gene")

# Save full results
write.csv(results, "results/tables/differential_expression_results.csv", row.names = FALSE)

# Save significant results
sig_results <- results %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1)

write.csv(sig_results, "results/tables/differential_expression_significant.csv", row.names = FALSE)

# Print summary
cat("Total genes tested:", nrow(results), "\n")
cat("Significant genes (adj.P.Val < 0.05 and |logFC| > 1):", nrow(sig_results), "\n")

# Show top genes
print(head(results, 10))