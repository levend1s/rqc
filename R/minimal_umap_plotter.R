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

all_cluster_levels <- levels(df$cluster)
cluster_to_show <- unique(df$cluster)

# Number of reads in the smallest sample
min_reads <- df %>%
  count(label) %>%
  summarise(min_n = min(n)) %>%
  pull(min_n)

# Downsample each sample to the same size
df_plot <- df %>%
  group_by(label) %>%
  slice_sample(n = min_reads) %>%
  ungroup()

df_plot %>%
  filter(cluster %in% cluster_to_show) %>%
  ggplot(aes(x = umap_x, y = umap_y, color = cluster)) +
  geom_point(size = 1, alpha = 0.6) +
  stat_density_2d(color = "lightgrey", linewidth = 0.5, ) +
  scale_color_discrete(drop = FALSE) +
  facet_wrap(~label) +
  labs(title = "UMAP by cluster, split by sample, downsampled for visual comparison between samples", x = "UMAP 1", y = "UMAP 2") +
  theme_classic()