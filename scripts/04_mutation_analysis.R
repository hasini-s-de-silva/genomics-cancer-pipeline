library(tidyverse)
library(ggplot2)

# -----------------------------
# Load mutation data
# -----------------------------
mutation_file <- "data/raw/mutation/tcga_luad_mutation.gz"

if (!file.exists(mutation_file)) {
  stop(paste("Mutation file not found at:", mutation_file))
}

cat("Reading mutation data...\n")

mutation_data <- read.delim(
  mutation_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

cat("Mutation data dimensions:", dim(mutation_data), "\n")

# -----------------------------
# Inspect structure
# -----------------------------
cat("Column names:\n")
print(colnames(mutation_data)[1:10])

# -----------------------------
# Convert to long format
# -----------------------------
# First column = genes
# Other columns = samples

mutation_long <- mutation_data %>%
  rename(gene = 1) %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "mutation_status"
  )

# Clean sample IDs
mutation_long$sample_id <- gsub("\\.", "-", mutation_long$sample_id)
mutation_long$sample_id <- gsub("^X", "", mutation_long$sample_id)
mutation_long$sample_id <- substr(mutation_long$sample_id, 1, 15)

# -----------------------------
# Filter mutated genes only
# -----------------------------
mutation_filtered <- mutation_long %>%
  filter(mutation_status != 0 & !is.na(mutation_status))

cat("Total mutation records:", nrow(mutation_filtered), "\n")

# -----------------------------
# Top mutated genes
# -----------------------------
top_mutated_genes <- mutation_filtered %>%
  group_by(gene) %>%
  summarise(mutation_count = n()) %>%
  arrange(desc(mutation_count)) %>%
  head(20)

print(top_mutated_genes)

# -----------------------------
# Plot: Top mutated genes
# -----------------------------
mutation_plot <- ggplot(top_mutated_genes,
                       aes(x = reorder(gene, mutation_count),
                           y = mutation_count)) +
  geom_bar(stat = "identity", fill = "#6A5ACD") +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    title = "Top Mutated Genes in TCGA LUAD",
    x = "Gene",
    y = "Mutation Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold")
  )

print(mutation_plot)

# Save plot
ggsave(
  "results/figures/top_mutated_genes.png",
  plot = mutation_plot,
  width = 7,
  height = 6
)

# -----------------------------
# Mutation burden per sample
# -----------------------------
mutation_burden <- mutation_filtered %>%
  group_by(sample_id) %>%
  summarise(total_mutations = n())

# Plot mutation burden
burden_plot <- ggplot(mutation_burden,
                      aes(x = total_mutations)) +
  geom_histogram(bins = 50, fill = "#00A087", color = "black") +
  theme_classic(base_size = 14) +
  labs(
    title = "Mutation Burden per Sample",
    x = "Number of Mutations",
    y = "Number of Samples"
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(burden_plot)

ggsave(
  "results/figures/mutation_burden.png",
  plot = burden_plot,
  width = 7,
  height = 5
)

cat("Mutation analysis complete.\n")