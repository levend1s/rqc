library(tidyverse)
library(umap)
library(randomForest)
library(dbscan)
library(ggrepel)
library(proxy)
library(processx)

all_continuous_cols <- c("read_start", "read_end", "read_length", "poly_a_length", "num_mods")
CLUSTER_PERC <- 0.01 # percentage of total reads to consider a cluster
UMAP_MIN_DIST <- 0.3

set.seed(42)

file <- c("~/rqc/cluster_transcripts_results.tsv")

# -------------------------------------------------------------------
# 8. Cluster summaries and labels
# -------------------------------------------------------------------

build_cluster_labels <- function(df,
                                 mod_matrix,
                                 continuous_cols,
                                 mod_pos_cols,
                                 top_n = 5) {
  
  overall_cont <- df %>%
    summarise(across(all_of(continuous_cols),
                     ~mean(.x, na.rm = TRUE))) %>%
    pivot_longer(
      everything(),
      names_to = "feature",
      values_to = "overall"
    )
  
  if (length(mod_pos_cols) > 0) {
    
    overall_mod <- mod_matrix %>%
      select(all_of(mod_pos_cols)) %>%
      summarise(across(everything(),
                       ~mean(.x, na.rm = TRUE))) %>%
      pivot_longer(
        everything(),
        names_to = "feature",
        values_to = "overall"
      )
    
  } else {
    
    overall_mod <- tibble(
      feature = character(),
      overall = numeric()
    )
    
  }
  
  overall <- bind_rows(overall_cont, overall_mod)
  
  clusters <- levels(df$cluster)
  clusters <- clusters[clusters != "0"]
  
  map_dfr(clusters, function(cl){
    
    ids <- df %>%
      filter(cluster == cl) %>%
      pull(read_id)
    
    cont_summary <- df %>%
      filter(cluster == cl) %>%
      summarise(across(all_of(continuous_cols),
                       ~mean(.x, na.rm = TRUE))) %>%
      pivot_longer(
        everything(),
        names_to = "feature",
        values_to = "value"
      )
    
    if(length(mod_pos_cols) > 0){
      
      mod_summary <- mod_matrix %>%
        filter(read_id %in% ids) %>%
        select(all_of(mod_pos_cols)) %>%
        summarise(across(everything(),
                         ~mean(.x, na.rm = TRUE))) %>%
        pivot_longer(
          everything(),
          names_to = "feature",
          values_to = "value"
        )
      
    } else {
      
      mod_summary <- tibble(
        feature = character(),
        value = numeric()
      )
      
    }
    
    summary <- bind_rows(cont_summary, mod_summary) %>%
      left_join(overall, by = "feature") %>%
      mutate(
        effect = abs(value - overall),
        direction = if_else(value > overall, "↑", "↓")
      ) %>%
      arrange(desc(effect))
    
    top_features <- head(summary, top_n)
    
    feature_text <- map_chr(
      seq_len(nrow(top_features)),
      function(i){
        
        f <- top_features$feature[i]
        v <- top_features$value[i]
        
        if(f %in% mod_pos_cols){
          
          sprintf("%s %.0f%%",
                  sub("^mod_pos_", "", f),
                  100*v)
          
        } else if(f %in% c("read_start",
                           "read_end",
                           "read_length")){
          
          sprintf("%s %.0f bp", f, v)
          
        } else if(f == "poly_a_length"){
          
          sprintf("%s %.1f nt", f, v)
          
        } else {
          
          sprintf("%s %.2f", f, v)
          
        }
        
      }
    )
    
    tibble(
      
      cluster = cl,
      
      n = length(ids),
      
      umap_x = median(df$umap_x[df$cluster==cl]),
      
      umap_y = median(df$umap_y[df$cluster==cl]),
      
      label = paste0(
        "Cluster ", cl,
        "\n(n=", length(ids), ")",
        "\n",
        paste(feature_text, collapse="\n")
      )
      
    )
    
  })
  
}

#---------------------------

df <- read.delim(file)

df <- df %>%
  mutate(
    mod_positions = str_remove_all(mod_positions, "\\[|\\]"),
    mod_positions = na_if(mod_positions, ""),
    mod_positions = str_split(mod_positions, ",\\s*")
  )


# -------------------------------------------------------------------
# 2. Build the mod_matrix (counts per position; presence/absence if a
#    position never repeats within a read's own list, which is the
#    typical case here)
# -------------------------------------------------------------------
df_long <- df %>%
  select(read_id, mod_positions) %>%
  unnest(mod_positions) %>%
  filter(!is.na(mod_positions),
         mod_positions != "") %>%
  rename(mod_position = mod_positions) %>%
  count(read_id, mod_position, name = "count")

mod_matrix <- df_long %>%
  pivot_wider(
    names_from = mod_position,
    values_from = count,
    values_fill = 0,
    names_prefix = "mod_pos_"
  )

mod_matrix <- df %>%
  select(read_id) %>%
  left_join(mod_matrix, by = "read_id") %>%
  replace(is.na(.), 0)

mod_matrix <- mod_matrix %>%
  mutate(
    no_mod_reads = as.numeric(
      rowSums(across(starts_with("mod_pos_"))) == 0
    )
  )

# DONT INCLUDE NO_MOD_READS
all_mod_pos_cols <- c(
  grep("^mod_pos_", names(mod_matrix), value = TRUE)
)
mod_pos_cols <- all_mod_pos_cols

mod_X <- mod_matrix %>% select(all_of(mod_pos_cols)) %>% as.matrix()

# -------------------------------------------------------------------
# 3. TF x frequency-weighting of the modification-position matrix
#    TF = relative frequency of a position within a read's own list
#    FW = "frequency weight" -- deliberately the INVERSE of classic
#         IDF. Positions present across many/most reads are treated
#         as the informative, consistent signal here (e.g. core
#         modification sites), so they get UP-weighted rather than
#         down-weighted. FW = log(1 + doc_freq), so a position in
#         every read gets close to log(1 + N) (max weight), and a
#         position seen in only 1-2 reads gets close to log(2)
#         (near-minimum weight).
# -------------------------------------------------------------------
compute_freq_weighted <- function(count_matrix) {
  row_sums <- rowSums(count_matrix)
  row_sums[row_sums == 0] <- 1  # reads with no mods: avoid div by 0
  tf <- count_matrix / row_sums
  
  doc_freq <- colSums(count_matrix > 0)
  freq_weight <- log(1 + doc_freq)  # higher for MORE common positions
  
  sweep(tf, 2, freq_weight, `*`)
}

mod_X_tfidf <- compute_freq_weighted(mod_X)

# -------------------------------------------------------------------
# 4. Continuous features -> scale, then Euclidean distance
# -------------------------------------------------------------------
continuous_cols <- all_continuous_cols
cont_X <- df %>%
  select(all_of(continuous_cols)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.x)))) %>%
  as.matrix()

# -------------------------------------------------------------------
# 5. Distance matrices, normalized, combined
# -------------------------------------------------------------------
# cosine distance on frequency-weighted modification profiles
# (replaces the earlier binary/Jaccard distance on raw presence;
#  common positions are up-weighted, not down-weighted -- see step 3)
tfidf_cosine_dist <- as.matrix(proxy::dist(mod_X_tfidf, method = "cosine"))
tfidf_cosine_dist[is.na(tfidf_cosine_dist)] <- 1

euclidean_dist <- as.matrix(dist(cont_X, method = "euclidean"))

normalize01 <- function(m) (m - min(m)) / (max(m) - min(m))

tfidf_cosine_dist_norm <- normalize01(tfidf_cosine_dist)
euclidean_dist_norm <- normalize01(euclidean_dist)

w_mod <- 0.2
w_cont <- 0.8

combined_dist <- w_mod * tfidf_cosine_dist_norm + w_cont * euclidean_dist_norm
cluster_sample_count <- nrow(df) * CLUSTER_PERC

# -------------------------------------------------------------------
# 6. Run UMAP ONCE on the combined data -> shared coordinate space
# -------------------------------------------------------------------
custom_config <- umap.defaults
custom_config$input <- "dist"
custom_config$n_neighbors <- cluster_sample_count
custom_config$min_dist <- UMAP_MIN_DIST
umap_result <- umap(combined_dist, config = custom_config)

df$umap_x <- umap_result$layout[, 1]
df$umap_y <- umap_result$layout[, 2]

# -------------------------------------------------------------------
# 7. Cluster ONCE on the combined distance matrix -> shared cluster IDs
# -------------------------------------------------------------------
clust <- hdbscan(umap_result$layout, minPts = cluster_sample_count)

df$cluster <- as.factor(clust$cluster)
table(df$cluster, df$source_file)  # check cluster x file distribution -- see batch-effect note below

# -------------------------------------------------------------------
# 9. Plot: same colors, same layout, same clusters, faceted by file
# -------------------------------------------------------------------
cluster_labels <- build_cluster_labels(
  df,
  mod_matrix,
  continuous_cols,
  mod_pos_cols,
  top_n = 5
)

all_cluster_levels <- levels(df$cluster)
cluster_to_show <- unique(df$cluster[df$cluster != "0"])
# cluster_to_show <- c("3", "4", "5", "6")
# cluster_to_show <- unique(df$cluster)

# Number of reads in the smallest sample
min_reads <- df %>%
  count(source_file) %>%
  summarise(min_n = min(n)) %>%
  pull(min_n)

# Downsample each sample to the same size
df_plot <- df %>%
  group_by(source_file) %>%
  slice_sample(n = min_reads) %>%
  ungroup()

df_plot %>%
  filter(cluster %in% cluster_to_show) %>%
  ggplot(aes(x = umap_x, y = umap_y, color = cluster)) +
  geom_point(size = 1, alpha = 0.6) +
  stat_density_2d(color = "lightgrey", linewidth = 0.5, ) +
  scale_color_discrete(drop = FALSE) +
  facet_wrap(~source_file) +
  labs(title = "UMAP by cluster, split by sample, downsampled for visual comparison between samples", x = "UMAP 1", y = "UMAP 2") +
  theme_classic()

# optional: single combined (non-faceted) plot with labels, for the
# overall cluster structure across all files at once
df %>%
  filter(cluster %in% cluster_to_show) %>%
  ggplot(aes(x = umap_x, y = umap_y, color = cluster)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_discrete(drop = FALSE) +
  geom_label_repel(
    data = cluster_labels %>% filter(cluster %in% cluster_to_show),
    aes(x = umap_x, y = umap_y, label = label),
    inherit.aes = FALSE,
    size = 2.8,
    fontface = "bold",
    fill = "white",
    alpha = 0.9,
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  labs(title = "UMAP colored by cluster (all samples combined)", x = "UMAP 1", y = "UMAP 2") +
  theme_classic() +
  theme(legend.position = "none")


