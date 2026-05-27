setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-MusMusculus-DRG-Pilot/pilot")

cat("\n==============================\n")
cat("Starting Volcano Plot Pipeline\n")
cat("==============================\n\n")

library(readr)
library(dplyr)
library(ggplot2)
library(tibble)

cat("Libraries loaded successfully.\n\n")

# =========================================================
# Create output directory
# =========================================================

output_dir <- "volcano_plot_analysis"

# if (!dir.exists(output_dir)) {
#   dir.create(output_dir)
# }

cat("Output directory ready:\n")
cat(output_dir, "\n\n")

# =========================================================
# 1. Load DESeq2 results
# =========================================================

cat("STEP 1: Loading DESeq2 results...\n")

res <- read_csv(
  "final_matrix/deseq2_results.csv"
)

cat("DESeq2 results loaded successfully.\n")
cat("Rows:", nrow(res), "\n")
cat("Columns:", ncol(res), "\n\n")

# =========================================================
# 2. Remove NA values
# =========================================================

cat("STEP 2: Removing NA values...\n")

res <- res %>%
  filter(
    !is.na(log2FoldChange),
    !is.na(padj)
  )

cat("Remaining genes after filtering:\n")
cat(nrow(res), "\n\n")

# =========================================================
# 3. Classify significance
# =========================================================

cat("STEP 3: Classifying significant genes...\n")

res <- res %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),

    neg_log10_padj = -log10(padj)
  )

cat("Significance classification complete.\n\n")

cat("Gene category counts:\n")
print(table(res$significance))
cat("\n")

# =========================================================
# 4. Save processed volcano results
# =========================================================

cat("STEP 4: Saving processed volcano results...\n")

write_csv(
  res,
  file.path(output_dir, "volcano_plot_results.csv")
)

cat("Saved processed volcano results.\n\n")

# =========================================================
# 5. Generate volcano plot
# =========================================================

cat("STEP 5: Generating volcano plot...\n")

volcano <- ggplot(
  res,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = significance
  )
) +
  geom_point(
    alpha = 0.7,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

cat("Volcano plot generated successfully.\n\n")

# =========================================================
# 6. Save volcano plot
# =========================================================

cat("STEP 6: Saving volcano plot...\n")

ggsave(
  file.path(output_dir, "volcano_plot.png"),
  volcano,
  width = 8,
  height = 6
)

cat("Saved volcano plot.\n\n")

# =========================================================
# 7. Generate summary statistics
# =========================================================

cat("STEP 7: Generating summary statistics...\n")

summary_table <- tibble(
  metric = c(
    "total_genes",
    "significant_genes",
    "upregulated_genes",
    "downregulated_genes"
  ),
  value = c(
    nrow(res),
    sum(res$padj < 0.05, na.rm = TRUE),
    sum(res$significance == "Upregulated"),
    sum(res$significance == "Downregulated")
  )
)

print(summary_table)

write_csv(
  summary_table,
  file.path(output_dir, "volcano_plot_summary.csv")
)

cat("Saved volcano plot summary.\n\n")

cat("==============================\n")
cat("VOLCANO PLOT PIPELINE COMPLETE\n")
cat("==============================\n")