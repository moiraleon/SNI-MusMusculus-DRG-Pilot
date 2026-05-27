setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-MusMusculus-DRG-Pilot/pilot")

cat("\n==============================\n")
cat("Starting PCA Plot Pipeline From Existing Matrix\n")
cat("==============================\n\n")

library(readr)
library(dplyr)
library(tibble)
library(ggplot2)

cat("Libraries loaded successfully.\n\n")

# =========================================================
# Create output directory
# =========================================================

output_dir <- "PCA_analysis"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

cat("Output directory ready:\n")
cat(output_dir, "\n\n")

# =========================================================
# 1. Load TPM matrix
# =========================================================

cat("STEP 1: Loading gene-level TPM matrix...\n")

tpm <- read_csv("final_matrix/gene_level_TPM_matrix.csv")

cat("TPM matrix loaded.\n")
cat("Rows:", nrow(tpm), "\n")
cat("Columns:", ncol(tpm), "\n\n")

# =========================================================
# 2. Load metadata
# =========================================================

cat("STEP 2: Loading sample metadata...\n")

metadata <- read_csv("metadata/sample_metadata.csv")

cat("Metadata loaded.\n")
cat("Samples:", nrow(metadata), "\n")
cat("Columns:\n")
print(colnames(metadata))
cat("\n")

# =========================================================
# 3. Prepare expression matrix
# =========================================================

cat("STEP 3: Preparing expression matrix...\n")

expr_matrix <- tpm %>%
  select(-gene_symbol) %>%
  column_to_rownames("ensembl_gene_id")

expr_matrix <- as.matrix(expr_matrix)
mode(expr_matrix) <- "numeric"

cat("Expression matrix prepared.\n")
cat("Genes:", nrow(expr_matrix), "\n")
cat("Samples:", ncol(expr_matrix), "\n\n")

# =========================================================
# 4. Match metadata and sample order
# =========================================================

cat("STEP 4: Matching metadata and sample order...\n")

metadata <- metadata %>%
  column_to_rownames("sample_id")

missing_in_expr <- setdiff(rownames(metadata), colnames(expr_matrix))
missing_in_metadata <- setdiff(colnames(expr_matrix), rownames(metadata))

if (length(missing_in_expr) > 0) {
  stop(
    "ERROR: These samples are in metadata but missing from TPM matrix:\n",
    paste(missing_in_expr, collapse = "\n")
  )
}

if (length(missing_in_metadata) > 0) {
  cat("WARNING: These samples are in TPM matrix but not metadata:\n")
  print(missing_in_metadata)
  cat("They will be removed from PCA.\n\n")
}

expr_matrix <- expr_matrix[, rownames(metadata)]

cat("Sample order matched successfully.\n\n")

# =========================================================
# 5. Log-transform TPM values
# =========================================================

cat("STEP 5: Log-transforming TPM values...\n")

log_expr <- log2(expr_matrix + 1)

cat("Log transformation complete.\n")
cat("Matrix dimensions:", dim(log_expr), "\n\n")

# =========================================================
# 6. Filter low-expression genes
# =========================================================

cat("STEP 6: Filtering low-expression genes...\n")

keep <- rowSums(log_expr > 1) >= 2
log_expr_filtered <- log_expr[keep, ]

cat("Genes retained:", nrow(log_expr_filtered), "\n")
cat("Genes removed:", sum(!keep), "\n\n")

# =========================================================
# 7. Run PCA
# =========================================================

cat("STEP 7: Running PCA...\n")

pca <- prcomp(
  t(log_expr_filtered),
  scale. = TRUE
)

percent_var <- round(
  100 * (pca$sdev^2 / sum(pca$sdev^2))
)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample_id") %>%
  left_join(
    metadata %>% rownames_to_column("sample_id"),
    by = "sample_id"
  )

cat("PCA complete.\n")
cat("PC1 variance:", percent_var[1], "%\n")
cat("PC2 variance:", percent_var[2], "%\n\n")

# =========================================================
# 8. Save PCA coordinates
# =========================================================

cat("STEP 8: Saving PCA coordinates...\n")

write_csv(
  pca_df,
  "PCA_analysis/PCA_from_TPM_coordinates.csv"
)

cat("Saved PCA coordinates.\n\n")

# =========================================================
# 9. Generate PCA plot by GSE
# =========================================================

cat("STEP 9: Generating PCA plot colored by GSE...\n")

pca_plot_gse <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = gse)
) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA from TPM Matrix Colored by GSE",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  )

ggsave(
  "PCA_analysis/PCA_from_TPM_by_GSE.png",
  pca_plot_gse,
  width = 8,
  height = 6
)

cat("Saved PCA plot colored by GSE.\n\n")

# =========================================================
# 10. Generate PCA plot by condition
# =========================================================

cat("STEP 10: Generating PCA plot colored by condition...\n")

pca_plot_condition <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = condition)
) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA from TPM Matrix Colored by Condition",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  )

ggsave(
  "PCA_analysis/PCA_from_TPM_by_condition.png",
  pca_plot_condition,
  width = 8,
  height = 6
)

cat("Saved PCA plot colored by condition.\n\n")

cat("==============================\n")
cat("PCA PLOT PIPELINE COMPLETE\n")
cat("==============================\n")