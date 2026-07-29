#!/usr/bin/env Rscript

library(ggplot2)

# change this path to your TSV file
infile <- c("~/rqc/cluster_transcripts_results.tsv.umap")
df <- read.delim(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# make cluster discrete for coloring
df$cluster <- as.factor(df$cluster)

ggplot(df, aes(x = umap_x, y = umap_y, color = cluster)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "UMAP colored by cluster",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cluster"
  )