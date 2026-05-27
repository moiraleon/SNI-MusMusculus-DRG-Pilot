setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-MusMusculus-DRG-Pilot/pilot")

# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages(c(
#   "readr",
#   "dplyr",
#   "tibble"
# ))

# #Plotting Packages
# install.packages(c(
#   "ggplot2",
#   "pheatmap"
# ))

# BiocManager::install(c(
#   "EnhancedVolcano"
# ))

cat("\n==============================\n")
cat("Starting DESeq2 Differential Expression Pipeline\n")
cat("==============================\n\n")

library(readr)
library(dplyr)
library(tibble)
library(DESeq2)

cat("Libraries loaded successfully.\n\n")

# =========================================================
# 1. Load gene-level counts matrix
# =========================================================

cat("STEP 1: Loading gene-level counts matrix...\n")

counts <- read_csv("final_matrix/gene_level_counts_matrix.csv")

cat("Counts matrix loaded.\n")
cat("Rows:", nrow(counts), "\n")
cat("Columns:", ncol(counts), "\n\n")

# =========================================================
# 2. Load sample metadata
# =========================================================

cat("STEP 2: Loading sample metadata...\n")

metadata <- read_csv("metadata/sample_metadata.csv")

cat("Metadata loaded.\n")
cat("Number of samples:", nrow(metadata), "\n")
cat("Metadata columns:\n")
print(colnames(metadata))
cat("\n")

if (!"sample_id" %in% colnames(metadata)) {
  stop("ERROR: metadata/sample_metadata.csv must contain a sample_id column.")
}

if (!"condition" %in% colnames(metadata)) {
  stop("ERROR: metadata/sample_metadata.csv must contain a condition column.")
}

# =========================================================
# 3. Separate gene info from count values
# =========================================================

cat("STEP 3: Separating gene identifiers and count matrix...\n")

gene_info <- counts %>%
  select(ensembl_gene_id, gene_symbol)

count_matrix <- counts %>%
  select(-gene_symbol) %>%
  column_to_rownames("ensembl_gene_id")

cat("Gene info and count matrix prepared.\n")
cat("Genes:", nrow(count_matrix), "\n")
cat("Samples in count matrix:", ncol(count_matrix), "\n\n")

# =========================================================
# 4. Prepare metadata and match sample order
# =========================================================

cat("STEP 4: Preparing metadata and matching sample order...\n")

metadata <- metadata %>%
  column_to_rownames("sample_id")

missing_in_counts <- setdiff(rownames(metadata), colnames(count_matrix))
missing_in_metadata <- setdiff(colnames(count_matrix), rownames(metadata))

if (length(missing_in_counts) > 0) {
  stop(
    "ERROR: These samples are in metadata but missing from count matrix:\n",
    paste(missing_in_counts, collapse = "\n")
  )
}

if (length(missing_in_metadata) > 0) {
  cat("WARNING: These samples are in count matrix but not metadata:\n")
  print(missing_in_metadata)
  cat("They will be removed from the DESeq2 analysis.\n\n")
}

count_matrix <- count_matrix[, rownames(metadata)]

cat("Sample order matched successfully.\n")
cat("Final count matrix dimensions:", dim(count_matrix), "\n\n")

# =========================================================
# 5. Clean and round counts
# =========================================================

cat("STEP 5: Cleaning and rounding count matrix...\n")

count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "numeric"

if (any(is.na(count_matrix))) {
  stop("ERROR: Count matrix contains NA values.")
}

count_matrix <- round(count_matrix)

cat("Counts converted to numeric matrix and rounded.\n")
cat("Minimum count:", min(count_matrix), "\n")
cat("Maximum count:", max(count_matrix), "\n\n")

# =========================================================
# 6. Create DESeq2 dataset
# =========================================================

cat("STEP 6: Creating DESeq2 dataset...\n")

metadata$condition <- factor(metadata$condition)
metadata$gse <- factor(metadata$gse)

cat("Condition by GSE table:\n")
print(table(metadata$gse, metadata$condition))
cat("\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ gse
#   design = ~ gse + condition
)

cat("DESeq2 dataset created successfully.\n")
cat("Genes before filtering:", nrow(dds), "\n\n")

# =========================================================
# 7. Filter low-count genes
# =========================================================

cat("STEP 7: Filtering low-count genes...\n")

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("Genes retained after filtering:", nrow(dds), "\n")
cat("Genes removed:", sum(!keep), "\n\n")

# =========================================================
# 8. Run DESeq2
# =========================================================

cat("STEP 8: Running DESeq2 differential expression analysis...\n")

dds <- DESeq(dds)

cat("DESeq2 analysis complete.\n\n")

# =========================================================
# 9. Extract results
# =========================================================

cat("STEP 9: Extracting DESeq2 results...\n")

res <- results(dds)

cat("Results extracted.\n")
cat("Result columns:\n")
print(colnames(as.data.frame(res)))
cat("\n")

res_df <- as.data.frame(res) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info, by = "ensembl_gene_id") %>%
  relocate(gene_symbol, .after = ensembl_gene_id)

cat("Gene symbols added to results.\n")
cat("Total result rows:", nrow(res_df), "\n\n")

# =========================================================
# 10. Save results
# =========================================================

cat("STEP 10: Saving DESeq2 results...\n")

write_csv(
  res_df,
  "final_matrix/deseq2_results.csv"
)

cat("Saved DESeq2 results to final_matrix/deseq2_results.csv\n\n")

# =========================================================
# 11. Summary
# =========================================================

cat("STEP 11: Generating summary statistics...\n")

summary_table <- tibble(
  metric = c(
    "total_genes_tested",
    "significant_genes_padj_0.05",
    "upregulated_log2FC_gt_1",
    "downregulated_log2FC_lt_minus_1"
  ),
  value = c(
    nrow(res_df),
    sum(res_df$padj < 0.05, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)
  )
)

#Naming for output CSV of DESeq2
write_csv(
  summary_table,
  "final_matrix/deseq2_summary.csv"
)

cat("Saved DESeq2 summary to final_matrix/deseq2_summary.csv\n\n")

print(summary_table)

cat("\n==============================\n")
cat("DESeq2 PIPELINE COMPLETE\n")
cat("==============================\n")