# diff_abundance_clusters.R
# install.packages("BiocManager")
# BiocManager::install("edgeR")
library(edgeR)
library(ggplot2)
library(dplyr)
library(ggrepel)

excl_infile <- c("~/rqc/./cluster_transcripts_results_batch_m6A_28hpi.tsv.cluster_summary.tsv")  # EDIT to your actual filename


# 1) Load counts table
infile <- c("~/rqc/cluster_transcripts_results_batch_m6A_28hpi.tsv")
counts <- read.delim(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
counts <- counts[!grepl("noise|NA|empty|MIT|API", counts$ID_cluster), ]
rownames(counts) <- counts$ID_cluster
counts$ID_cluster <- NULL

# 2) Define sample groups (edit to your design)
group <- factor(c("C","C","K","K"))
stopifnot(length(group) == ncol(counts))

# 3) Build edgeR object
y <- DGEList(counts = counts, group = group)

# 4) Normalize + estimate dispersion
y <- calcNormFactors(y)
design <- model.matrix(~ group)
y <- estimateDisp(y, design)

# 5) Fit and test
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

# 6) Results
res <- topTags(qlf, n = Inf)$table
res$ID_cluster <- rownames(res)
res <- res[, c("ID_cluster", "logFC","logCPM","F","PValue","FDR")]
res <- res %>%
  mutate(
    ID      = sub("^([^_]*_[^_]*)_.*$", "\\1", ID_cluster),
    cluster = sub("^[^_]*_[^_]*_", "", ID_cluster)
  )
res$negLog10FDR <- -log10(res$FDR + 1e-300)
res$signif <- res$FDR < 0.05 & abs(res$logFC) > 1

# --- 7) Load exclusivity annotation file and merge on ID_cluster ---
excl <- read.delim(excl_infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
excl$ID_cluster <- excl$index

# combined_exclusivity_score may be stored as a stringified list (e.g. "[0.374]"
# or "[]") since it originates from a Python list column written to TSV.
# Treat empty/NA/empty-list values as "no score".
excl <- excl %>%
  mutate(
    has_exclusivity_score = !(is.na(combined_exclusivity_score) |
                                combined_exclusivity_score %in% c("", "[]", "NA"))
  ) %>%
  select(ID_cluster, has_exclusivity_score, target, exclusive_predictors)

res <- res %>%
  left_join(excl, by = "ID_cluster") %>%
  mutate(has_exclusivity_score = tidyr::replace_na(has_exclusivity_score, FALSE))

# --- 8) Combine significance + exclusivity into one coloring variable ---
res <- res %>%
  mutate(
    plot_group = case_when(
      signif & has_exclusivity_score  ~ "Significant + has exclusivity score",
      signif & !has_exclusivity_score ~ "Significant",
      !signif & has_exclusivity_score ~ "Has exclusivity score",
      TRUE                             ~ "Not significant"
    )
  )

top <- head(res[order(res$FDR), ], 15)

p <- ggplot(res, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = plot_group), alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "Not significant"                      = "grey70",
    "Significant"                          = "red3",
    "Has exclusivity score"                = "steelblue",
    "Significant + has exclusivity score"  = "purple3"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = top,
    aes(label = ID_cluster),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0
  ) +
  labs(
    title = "Differential abundance volcano plot",
    x = "log2 Fold Change (K vs C)",
    y = "-log10(FDR)",
    color = NULL
  ) +
  theme_bw() +
  theme(legend.position = "right")

print(p)

a <- res[res$has_exclusivity_score == TRUE,] %>% select(ID, target, exclusive_predictors)
