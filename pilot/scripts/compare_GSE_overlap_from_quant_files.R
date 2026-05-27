setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-MusMusculus-DRG-Pilot/pilot")

cat("\n==============================\n")
cat("Starting Expression Comparison Pipeline\n")
cat("==============================\n\n")

library(readr)
library(dplyr)
library(tibble)
library(tximport)
library(biomaRt)

cat("Libraries loaded successfully.\n\n")

# =========================================================
# 1. Load sample metadata
# =========================================================

cat("STEP 1: Loading sample metadata...\n")

samples <- read_csv("metadata/sample_metadata.csv")

cat("Metadata loaded.\n")
cat("Number of samples:", nrow(samples), "\n\n")

# =========================================================
# 2. Locate Salmon quant files
# =========================================================

cat("STEP 2: Locating Salmon quant.sf files...\n")

files <- file.path(
  "quant",
  samples$gse,
  samples$sample_id,
  "quant.sf"
)

names(files) <- samples$sample_id

cat("Constructed quant.sf paths:\n")
print(files)

cat("\nChecking file existence...\n")
print(file.exists(files))

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {
  stop(
    "\nERROR: Missing quant.sf files:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("All quant.sf files located successfully.\n\n")

# =========================================================
# 3. Build transcript-to-gene mapping
# =========================================================

cat("STEP 3: Connecting to Ensembl BioMart mirror...\n")

mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror = "asia"
)

cat("Connected to Ensembl.\n")
cat("Downloading transcript-to-gene mapping...\n")

tx2gene <- getBM(
  attributes = c(
    "ensembl_transcript_id_version",
    "ensembl_gene_id",
    "external_gene_name"
  ),
  mart = mart
)

cat("Transcript-to-gene mapping downloaded.\n")
cat("Rows in tx2gene:", nrow(tx2gene), "\n")

write_csv(
  tx2gene,
  "gene_mapping/ensembl_tx2gene_mapping.csv"
)

cat("Saved transcript-to-gene mapping.\n\n")

tx2gene_simple <- tx2gene %>%
  dplyr::select(
    ensembl_transcript_id_version,
    ensembl_gene_id
  ) %>%
  distinct()

cat("Simplified tx2gene mapping created.\n\n")

# =========================================================
# 4. Import Salmon results
# =========================================================

cat("STEP 4: Importing Salmon quantification files...\n")

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene_simple,
  ignoreTxVersion = FALSE
)

cat("Salmon quantifications imported successfully.\n")
cat("Genes imported:", nrow(txi$abundance), "\n")
cat("Samples imported:", ncol(txi$abundance), "\n\n")

# =========================================================
# 5. Create TPM and counts matrices
# =========================================================

cat("STEP 5: Creating TPM and counts matrices...\n")

gene_tpm <- as.data.frame(txi$abundance) %>%
  rownames_to_column("ensembl_gene_id")

gene_counts <- as.data.frame(txi$counts) %>%
  rownames_to_column("ensembl_gene_id")

cat("TPM matrix dimensions:", dim(gene_tpm), "\n")
cat("Counts matrix dimensions:", dim(gene_counts), "\n\n")

# =========================================================
# 6. Add gene symbols
# =========================================================

cat("STEP 6: Adding gene symbols...\n")

gene_symbols <- tx2gene %>%
  dplyr::select(
    ensembl_gene_id,
    external_gene_name
  ) %>%
  distinct() %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    gene_symbol = first(na.omit(external_gene_name)),
    .groups = "drop"
  )

gene_tpm <- gene_tpm %>%
  left_join(gene_symbols, by = "ensembl_gene_id") %>%
  relocate(gene_symbol, .after = ensembl_gene_id)

gene_counts <- gene_counts %>%
  left_join(gene_symbols, by = "ensembl_gene_id") %>%
  relocate(gene_symbol, .after = ensembl_gene_id)

cat("Gene symbols added successfully.\n\n")

# =========================================================
# 7. Save TPM and counts matrices
# =========================================================

cat("STEP 7: Saving gene-level matrices...\n")

write_csv(
  gene_tpm,
  "final_matrix/gene_level_TPM_matrix.csv"
)

write_csv(
  gene_counts,
  "final_matrix/gene_level_counts_matrix.csv"
)

cat("Saved TPM matrix.\n")
cat("Saved counts matrix.\n\n")

# =========================================================
# 8. Expression overlap analysis
# =========================================================

cat("STEP 8: Performing overlap analysis...\n")

tpm_matrix <- gene_tpm %>%
  dplyr::select(-gene_symbol) %>%
  column_to_rownames("ensembl_gene_id")

expressed <- tpm_matrix > 1

cat("Expression threshold set to TPM > 1\n")

gse261676_samples <- samples %>%
  filter(gse == "GSE261676") %>%
  pull(sample_id)

gse123919_samples <- samples %>%
  filter(gse == "GSE123919") %>%
  pull(sample_id)

cat("GSE261676 samples:", length(gse261676_samples), "\n")
cat("GSE123919 samples:", length(gse123919_samples), "\n")

overlap_table <- tibble(
  ensembl_gene_id = rownames(tpm_matrix),

  expressed_in_GSE261676 =
    rowSums(
      expressed[, gse261676_samples, drop = FALSE]
    ) > 0,

  expressed_in_GSE123919 =
    rowSums(
      expressed[, gse123919_samples, drop = FALSE]
    ) > 0,

  n_samples_expressed_GSE261676 =
    rowSums(
      expressed[, gse261676_samples, drop = FALSE]
    ),

  n_samples_expressed_GSE123919 =
    rowSums(
      expressed[, gse123919_samples, drop = FALSE]
    )

) %>%
  left_join(
    gene_symbols,
    by = "ensembl_gene_id"
  ) %>%
  relocate(
    gene_symbol,
    .after = ensembl_gene_id
  ) %>%
  mutate(
    overlap_status = case_when(
      expressed_in_GSE261676 &
        expressed_in_GSE123919 ~ "shared",

      expressed_in_GSE261676 &
        !expressed_in_GSE123919 ~ "GSE261676_only",

      !expressed_in_GSE261676 &
        expressed_in_GSE123919 ~ "GSE123919_only",

      TRUE ~ "not_expressed"
    )
  )

cat("Overlap analysis complete.\n\n")

write_csv(
  overlap_table,
  "final_matrix/expression_overlap_table.csv"
)

cat("Saved overlap table.\n")

shared_genes <- overlap_table %>%
  filter(overlap_status == "shared")

write_csv(
  shared_genes,
  "final_matrix/shared_expressed_genes.csv"
)

cat("Saved shared expressed genes table.\n\n")

# =========================================================
# 9. Summary
# =========================================================

cat("STEP 9: Generating summary statistics...\n")

summary_table <- overlap_table %>%
  count(overlap_status)

write_csv(
  summary_table,
  "final_matrix/expression_overlap_summary.csv"
)

cat("Saved overlap summary.\n\n")

cat("==============================\n")
cat("PIPELINE COMPLETE\n")
cat("==============================\n\n")

print(summary_table)