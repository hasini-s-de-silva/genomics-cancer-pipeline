library(tidyverse)

# Load expression data
expr <- read.delim("data/raw/expression/tcga_luad_expression.tsv", check.names = FALSE)

# Extract sample IDs (column names)
samples <- colnames(expr)[-1]

# Create metadata dataframe
metadata <- data.frame(
  sample_id = samples
)

# Extract sample type from barcode
metadata$sample_type_code <- substr(metadata$sample_id, 14, 15)

# Map to labels
metadata$sample_type <- ifelse(metadata$sample_type_code == "01", "Tumor",
                        ifelse(metadata$sample_type_code == "11", "Normal", "Other"))

# Check counts
table(metadata$sample_type)

# Save metadata
write.csv(metadata, "data/processed/metadata_clean.csv", row.names = FALSE)

# Preview
head(metadata)