library(tidyverse)
library(umap)
library(randomForest)
library(dbscan)
library(ggrepel)
library(proxy)


# SAMPLE=28C1; python rqc.py cluster_transcripts -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --padding 100 --ids PF3D7_1442300.1 -b /Users/joshualevendis/Downloads/bams/${SAMPLE}_to_pfal.50MAPQ.sorted.bam -o ~/Downloads/umap_results_${SAMPLE}.tsv

# =====================================================================
# CONFIG
# =====================================================================

all_continuous_cols <- c("read_start", "read_end", "read_length", "poly_a_length", "num_mods")
MIN_PTS_PCT <- 0.05
MINPTS_CLUSTER <- 10

set.seed(42)


file <- c("~/rqc/cluster_transcripts_results.tsv")

df <- read.delim(file)

df <- df %>%
  mutate(
    mod_positions = str_remove_all(mod_positions, "\\[|\\]"),
    mod_positions = na_if(mod_positions, ""),
    mod_positions = str_split(mod_positions, ",\\s*")
  )

# -------------------------------------------------------------------
# 2. Build the binary mod_matrix (presence/absence per position)
# -------------------------------------------------------------------
df_long <- df %>%
  select(read_id, mod_positions) %>%
  unnest(mod_positions) %>%
  filter(!is.na(mod_positions),
         mod_positions != "") %>%
  rename(mod_position = mod_positions) %>%
  mutate(present = 1)

mod_matrix <- df_long %>%
  select(read_id, mod_position, present) %>%
  pivot_wider(
    names_from = mod_position,
    values_from = present,
    values_fill = 0,
    names_prefix = "mod_pos_"
  )

mod_matrix <- df %>%
  select(read_id) %>%
  left_join(mod_matrix, by="read_id") %>%
  replace(is.na(.), 0)

mod_matrix <- mod_matrix %>%
  mutate(
    no_mod_reads = as.numeric(
      rowSums(across(starts_with("mod_pos_"))) == 0
    )
  )
# all_mod_pos_cols <- c(
#   grep("^mod_pos_", names(mod_matrix), value = TRUE),
#   "no_mod_reads"
# )

# DONT INCLUDE NO_MOD_READS
all_mod_pos_cols <- c(
  grep("^mod_pos_", names(mod_matrix), value = TRUE)
)

mod_pos_cols <- all_mod_pos_cols

mod_X <- mod_matrix %>% select(all_of(mod_pos_cols)) %>% as.matrix()

# -------------------------------------------------------------------
# 3. Continuous features -> scale, then Euclidean distance
# -------------------------------------------------------------------
continuous_cols <- all_continuous_cols

cont_X <- df %>%
  select(all_of(continuous_cols)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.x)))) %>%
  as.matrix()

# -------------------------------------------------------------------
# 4. Distance matrices, normalized, combined
# -------------------------------------------------------------------
# jaccard_dist <- as.matrix(dist(mod_X, method = "binary"))
# jaccard_dist[is.na(jaccard_dist)] <- 1

# cosine distance for modification profiles
cosine_dist <- as.matrix(proxy::dist(mod_X, method = "cosine"))
cosine_dist[is.na(cosine_dist)] <- 1

euclidean_dist <- as.matrix(dist(cont_X, method = "euclidean"))

normalize01 <- function(m) (m - min(m)) / (max(m) - min(m))
# jaccard_dist_norm <- normalize01(jaccard_dist)
cosine_dist_norm <- normalize01(cosine_dist)
euclidean_dist_norm <- normalize01(euclidean_dist)

w_mod <- 0.2
w_cont <- 0.8
# combined_dist <- w_mod * jaccard_dist_norm + w_cont * euclidean_dist_norm
combined_dist <- w_mod * cosine_dist_norm + w_cont * euclidean_dist_norm

# -------------------------------------------------------------------
# 5. Run UMAP ONCE on the combined data -> shared coordinate space
# -------------------------------------------------------------------
custom_config <- umap.defaults
custom_config$input <- "dist"
# custom_config$n_neighbors <- 30
# custom_config$min_dist <- 0.3
# custom_config$spread <- 1.5
umap_result <- umap(combined_dist, config = custom_config)

df$umap_x <- umap_result$layout[, 1]
df$umap_y <- umap_result$layout[, 2]

# -------------------------------------------------------------------
# 6. Cluster ONCE on the combined distance matrix -> shared cluster IDs
# -------------------------------------------------------------------
min_pts <- max(2, round(nrow(df) / length(unique(df$source_file)) * MIN_PTS_PCT))
clust <- hdbscan(as.dist(combined_dist), minPts = MINPTS_CLUSTER)
df$cluster <- as.factor(clust$cluster)

table(df$cluster, df$source_file)  # check cluster x file distribution -- see batch-effect note below


# -------------------------------------------------------------------
# 9. Cluster labels (top 5 distinguishing features per cluster)
# -------------------------------------------------------------------

build_cluster_labels <- function(df, mod_matrix,
                                 continuous_cols = c("read_end",
                                                     "read_length",
                                                     "poly_a_length",
                                                     "read_start"),
                                 mod_pos_cols,
                                 top_n_features = 5) {
  
  clusters <- df %>%
    filter(cluster != "0") %>%
    pull(cluster) %>%
    unique() %>%
    sort()
  
  # Overall feature means
  overall_mod <- mod_matrix %>%
    select(all_of(mod_pos_cols)) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(
      everything(),
      names_to = "feature",
      values_to = "overall_value"
    )
  
  overall_cont <- df %>%
    summarise(across(all_of(continuous_cols),
                     ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(
      everything(),
      names_to = "feature",
      values_to = "overall_value"
    )
  
  overall_all <- bind_rows(overall_mod, overall_cont)
  
  
  # Format labels
  format_line <- function(feature, value, is_mod) {
    
    if (is_mod) {
      
      feature_name <- sub("^mod_pos_", "", feature)
      
      return(
        paste0(feature_name,
               ": ",
               round(value * 100),
               "%")
      )
      
    } else if (feature == "read_end") {
      
      return(paste0(feature, ": ", round(value), " bp"))
      
    } else if (feature == "read_start") {
      
      return(paste0(feature, ": ", round(value), " bp"))
      
    } else if (feature == "read_length") {
      
      return(paste0(feature, ": ", round(value), " bp"))
      
    } else if (feature == "poly_a_length") {
      
      return(paste0(feature, ": ", round(value, 1), " nt"))
      
    } else {
      
      return(paste0(feature, ": ", round(value, 2)))
    }
  }
  
  
  map_dfr(clusters, function(cl) {
    
    ids <- df %>%
      filter(cluster == cl) %>%
      pull(read_id)
    
    
    centroid_x <- median(df$umap_x[df$cluster == cl])
    centroid_y <- median(df$umap_y[df$cluster == cl])
    
    
    # Cluster feature values
    this_mod <- mod_matrix %>%
      filter(read_id %in% ids) %>%
      select(all_of(mod_pos_cols)) %>%
      summarise(across(everything(), mean)) %>%
      pivot_longer(
        everything(),
        names_to = "feature",
        values_to = "cluster_value"
      )
    
    
    this_cont <- df %>%
      filter(cluster == cl) %>%
      summarise(across(all_of(continuous_cols),
                       ~ mean(.x, na.rm = TRUE))) %>%
      pivot_longer(
        everything(),
        names_to = "feature",
        values_to = "cluster_value"
      )
    
    
    # Combine all features
    this_all <- bind_rows(this_mod, this_cont) %>%
      left_join(overall_all, by = "feature") %>%
      mutate(
        is_mod = feature %in% mod_pos_cols,
        
        # absolute difference from global average
        effect = abs(cluster_value - overall_value)
      )
    
    
    # Select top 5 features describing this cluster
    top_features <- this_all %>%
      arrange(desc(effect)) %>%
      slice_head(n = top_n_features)
    
    
    feature_lines <- top_features %>%
      mutate(
        line = mapply(
          format_line,
          feature,
          cluster_value,
          is_mod
        )
      ) %>%
      pull(line)
    
    
    label_text <- paste0(
      "Cluster ", cl,
      " (n=", length(ids), ")\n",
      paste(feature_lines, collapse = "\n")
    )
    
    
    tibble(
      cluster = cl,
      umap_x = centroid_x,
      umap_y = centroid_y,
      label = label_text
    )
    
  })
}


cluster_labels <- build_cluster_labels(
  df,
  mod_matrix,
  continuous_cols = c("read_end",
                      "read_length",
                      "poly_a_length",
                      "read_start"),
  mod_pos_cols,
  top_n_features = 5
)
# -------------------------------------------------------------------
# 10. Plot: same colors, same layout, same clusters, faceted by file
# -------------------------------------------------------------------
all_cluster_levels <- levels(df$cluster)

cluster_to_show <- unique(df$cluster[df$cluster != "0"])
# cluster_to_show <- c("3", "4", "5", "6")
cluster_to_show <- unique(df$cluster)

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

df %>%
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

# -------------------------------------------------------------------
# 8. Random forest importance (on combined data)
# -------------------------------------------------------------------
rf_data <- df %>%
  select(read_id, cluster, all_of(continuous_cols)) %>%
  left_join(mod_matrix %>% select(read_id, all_of(mod_pos_cols)), by = "read_id") %>%
  select(-read_id) %>%
  drop_na()

rf_data <- rf_data %>% filter(cluster != "0")
rf_data$cluster <- droplevels(rf_data$cluster)

rf <- randomForest(cluster ~ ., data = rf_data, importance = TRUE)
importance(rf) %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  head(20)
