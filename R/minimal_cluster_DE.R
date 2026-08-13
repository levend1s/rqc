# diff_abundance_clusters.R
# install.packages("BiocManager")
# BiocManager::install("edgeR")

library(edgeR)
library(ggplot2)


# 1) Load counts table
# File format: tab-delimited with first column = ID_cluster
infile <- c("~/rqc/cluster_transcripts_results.tsv")
infile <- c("~/rqc/cluster_transcripts_results_batch_m6A.tsv")

counts <- read.delim(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
counts <- counts[!grepl("noise|NA|empty|MIT|API", counts$ID_cluster), ]


rownames(counts) <- counts$ID_cluster
counts$ID_cluster <- NULL

# 2) Define sample groups (edit to your design)
# Here: 28C* = C, 28K* = K
group <- factor(c("C","C","K","K"))
stopifnot(length(group) == ncol(counts))

# 3) Build edgeR object
y <- DGEList(counts = counts, group = group)

# Optional: filter very low-abundance rows
# keep <- filterByExpr(y, group = group)
# y <- y[keep, , keep.lib.sizes = FALSE]

# 4) Normalize + estimate dispersion
y <- calcNormFactors(y)
design <- model.matrix(~ group)   # tests K vs C
y <- estimateDisp(y, design)

# 5) Fit and test
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)  # groupK effect

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

top <- head(res[order(res$FDR), ], 15)

p <- ggplot(res, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = signif), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red3")) +
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
    y = "-log10(FDR)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# ggsave("cluster_DA_volcano.png", p, width = 6, height = 5, dpi = 300)
print(p)

