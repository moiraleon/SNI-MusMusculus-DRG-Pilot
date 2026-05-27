setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-MusMusculus-DRG-Pilot/pilot")

#Generate form batch corrected TPM matrix - using : "compare_GSE_overlap_from_quant_files.R"
cat("\n==============================\n")
cat("Starting Heatmap Generation Pipeline\n")
cat("==============================\n\n")

install.packages("pheatmap")

library(readr)
library(dplyr)
library(tibble)
library(pheatmap)
library(limma)

cat("Libraries loaded successfully.\n\n")

# =========================================================
# Create output directory
# =========================================================

output_dir <- "heatmap_analysis"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

cat("Output directory ready:", output_dir, "\n\n")

# =========================================================
# 1. Load TPM matrix
# =========================================================

cat("STEP 1: Loading TPM matrix...\n")

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
cat("Samples:", nrow(metadata), "\n\n")

# =========================================================
# 3. Prepare expression matrix
# =========================================================

cat("STEP 3: Preparing expression matrix...\n")

gene_symbols <- tpm %>%
  select(ensembl_gene_id, gene_symbol)

expr_matrix <- tpm %>%
  select(-gene_symbol) %>%
  column_to_rownames("ensembl_gene_id")

expr_matrix <- as.matrix(expr_matrix)
mode(expr_matrix) <- "numeric"

metadata <- metadata %>%
  column_to_rownames("sample_id")

expr_matrix <- expr_matrix[, rownames(metadata)]

cat("Expression matrix prepared.\n")
cat("Genes:", nrow(expr_matrix), "\n")
cat("Samples:", ncol(expr_matrix), "\n\n")

# =========================================================
# 4. Log-transform TPM
# =========================================================

cat("STEP 4: Log-transforming TPM values...\n")

log_expr <- log2(expr_matrix + 1)

cat("Log transformation complete.\n\n")

# =========================================================
# 5. Remove batch effect
# =========================================================

cat("STEP 5: Removing batch effect using GSE...\n")

batch_corrected <- removeBatchEffect(
  log_expr,
  batch = metadata$gse
)

cat("Batch correction complete.\n\n")

# =========================================================
# 6. Select top variable genes
# =========================================================

cat("STEP 6: Selecting top variable genes...\n")

gene_variance <- apply(batch_corrected, 1, var)

top_genes <- names(
  sort(gene_variance, decreasing = TRUE)
)[1:500]

heatmap_matrix <- batch_corrected[top_genes, ]

cat("Top variable genes selected:", length(top_genes), "\n\n")

# =========================================================
# 7. Scale genes
# =========================================================

cat("STEP 7: Scaling genes...\n")

heatmap_scaled <- t(scale(t(heatmap_matrix)))

cat("Scaling complete.\n\n")

# =========================================================
# 8. Prepare annotation
# =========================================================

cat("STEP 8: Preparing sample annotation...\n")

annotation_col <- metadata %>%
  select(gse, condition)

cat("Annotation prepared.\n\n")

# =========================================================
# 9. Generate heatmap
# =========================================================

cat("STEP 9: Generating heatmap...\n")

png(
  file.path(output_dir, "heatmap_top_500_variable_genes_batch_corrected.png"),
  width = 1200,
  height = 1000
)

pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 500 Variable Genes After Batch Correction"
)

dev.off()

cat("Heatmap saved.\n\n")

# =========================================================
# 10. Save heatmap matrix
# =========================================================

cat("STEP 10: Saving heatmap matrix...\n")

heatmap_df <- as.data.frame(heatmap_scaled) %>%
  rownames_to_column("ensembl_gene_id")

write_csv(
  heatmap_df,
  file.path(output_dir, "heatmap_top_500_variable_genes_matrix.csv")
)

cat("Heatmap matrix saved.\n\n")

cat("==============================\n")
cat("HEATMAP PIPELINE COMPLETE\n")
cat("==============================\n")