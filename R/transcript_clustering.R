library(tidyverse)
library(umap)
library(randomForest)
library(dbscan)
library(ggrepel)


# SAMPLE=28C1; python rqc.py cluster_transcripts -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --padding 100 --ids PF3D7_1442300.1 -b /Users/joshualevendis/Downloads/bams/${SAMPLE}_to_pfal.50MAPQ.sorted.bam -o ~/Downloads/umap_results_${SAMPLE}.tsv

# =====================================================================
# CONFIG
# =====================================================================
INCLUDE_COLS <- c("read_end", "read_length", "read_start", "poly_a_length", "read_strand")
INCLUDE_MOD_POS <- c("mod_pos_946343", "mod_pos_946455", "mod_pos_946461")
INCLUDE_MOD_POS <- NULL
MIN_PTS_PCT <- 0.05

set.seed(42)


files <- c(
  "~/Downloads/umap_results_28C1.tsv",
  "~/Downloads/umap_results_28C2.tsv",
  "~/Downloads/umap_results_28K1.tsv",
  "~/Downloads/umap_results_28K2.tsv"
)

# =====================================================================
# 1. Load and combine all files, tagging each row's source
# =====================================================================
df <- map_dfr(files, function(f) {
  read_tsv(f) %>% mutate(source_file = str_remove(basename(f), "\\.tsv$"))
})

df <- df %>%
  mutate(mod_positions = str_remove_all(mod_positions, "\\[|\\]")) %>%
  mutate(mod_positions = map(mod_positions, ~ as.numeric(str_split(.x, ",\\s*")[[1]])))

df <- df %>% mutate(row_id = row_number())

# -------------------------------------------------------------------
# 2. Build the binary mod_matrix (presence/absence per position)
# -------------------------------------------------------------------
df_long <- df %>%
  unnest(mod_positions) %>%
  rename(mod_position = mod_positions) %>%
  mutate(present = 1)

mod_matrix <- df_long %>%
  select(row_id, mod_position, present) %>%
  pivot_wider(
    names_from = mod_position,
    values_from = present,
    values_fill = 0,
    names_prefix = "mod_pos_"
  )

mod_matrix <- df %>%
  select(row_id) %>%
  left_join(mod_matrix, by = "row_id") %>%
  replace(is.na(.), 0)

all_mod_pos_cols <- grep("^mod_pos_", names(mod_matrix), value = TRUE)

mod_pos_cols <- if (is.null(INCLUDE_MOD_POS)) {
  all_mod_pos_cols
} else {
  intersect(all_mod_pos_cols, INCLUDE_MOD_POS)
}

mod_X <- mod_matrix %>% select(all_of(mod_pos_cols)) %>% as.matrix()

# -------------------------------------------------------------------
# 3. Continuous features -> scale, then Euclidean distance
# -------------------------------------------------------------------
all_continuous_cols <- c("read_start", "read_end", "read_length", "poly_a_length")
continuous_cols <- intersect(all_continuous_cols, INCLUDE_COLS)

cont_X <- df %>%
  select(all_of(continuous_cols)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.x)))) %>%
  as.matrix()

# -------------------------------------------------------------------
# 4. Distance matrices, normalized, combined
# -------------------------------------------------------------------
jaccard_dist <- as.matrix(dist(mod_X, method = "binary"))
euclidean_dist <- as.matrix(dist(cont_X, method = "euclidean"))

jaccard_dist[is.na(jaccard_dist)] <- 1

normalize01 <- function(m) (m - min(m)) / (max(m) - min(m))
jaccard_dist_norm <- normalize01(jaccard_dist)
euclidean_dist_norm <- normalize01(euclidean_dist)

w_mod <- 0.5
w_cont <- 0.5
combined_dist <- w_mod * jaccard_dist_norm + w_cont * euclidean_dist_norm

# -------------------------------------------------------------------
# 5. Run UMAP ONCE on the combined data -> shared coordinate space
# -------------------------------------------------------------------
custom_config <- umap.defaults
custom_config$input <- "dist"
custom_config$n_neighbors <- 30
custom_config$min_dist <- 0.3
custom_config$spread <- 1.5
umap_result <- umap(combined_dist, config = custom_config)

df$umap_x <- umap_result$layout[, 1]
df$umap_y <- umap_result$layout[, 2]

# -------------------------------------------------------------------
# 6. Cluster ONCE on the combined distance matrix -> shared cluster IDs
# -------------------------------------------------------------------
min_pts <- max(2, round(nrow(df) / length(files) * MIN_PTS_PCT))
clust <- hdbscan(as.dist(combined_dist), minPts = min_pts)
df$cluster <- as.factor(clust$cluster)

table(df$cluster, df$source_file)  # check cluster x file distribution -- see batch-effect note below

# -------------------------------------------------------------------
# 7. Batch-effect diagnostic: is source_file predicting cluster too well?
# -------------------------------------------------------------------
batch_check <- df %>%
  filter(cluster != "0") %>%
  count(cluster, source_file) %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n) * 100) %>%
  arrange(cluster, desc(pct))

print(batch_check, n = Inf)
# if any cluster is dominated (e.g. >90%) by a single source_file,
# that cluster may reflect a batch effect rather than real biology --
# worth extra scrutiny before drawing conclusions from it

# -------------------------------------------------------------------
# 8. Random forest importance (on combined data)
# -------------------------------------------------------------------
# rf_data <- df %>%
#   select(row_id, cluster, all_of(continuous_cols)) %>%
#   left_join(mod_matrix %>% select(row_id, all_of(mod_pos_cols)), by = "row_id") %>%
#   select(-row_id) %>%
#   drop_na()
# 
# rf_data <- rf_data %>% filter(cluster != "0")
# rf_data$cluster <- droplevels(rf_data$cluster)
# 
# rf <- randomForest(cluster ~ ., data = rf_data, importance = TRUE)
# importance(rf) %>%
#   as.data.frame() %>%
#   rownames_to_column("feature") %>%
#   arrange(desc(MeanDecreaseGini)) %>%
#   head(20)

# -------------------------------------------------------------------
# 9. Cluster labels (computed once, shared across all facets)
# -------------------------------------------------------------------
build_cluster_labels <- function(df, mod_matrix,
                                 continuous_cols = c("read_end", "read_length", "poly_a_length", "read_start"),
                                 mod_pos_cols, top_n_mods = 2) {
  clusters <- df %>% filter(cluster != "0") %>% pull(cluster) %>% unique() %>% sort()
  
  overall_mod_prev <- mod_matrix %>%
    select(all_of(mod_pos_cols)) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(everything(), names_to = "feature", values_to = "overall_value")
  
  overall_cont_mean <- df %>%
    summarise(across(all_of(continuous_cols), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "feature", values_to = "overall_value")
  
  overall_all <- bind_rows(overall_mod_prev, overall_cont_mean)
  
  format_line <- function(feature, value, is_mod) {
    case_when(
      is_mod ~ paste0(feature, ": ", round(value * 100), "%"),
      feature == "read_end" ~ paste0(feature, ": ", round(value), " bp"),
      feature == "read_start" ~ paste0(feature, ": ", round(value), " bp"),
      feature == "read_length" ~ paste0(feature, ": ", round(value), " bp"),
      feature == "poly_a_length" ~ paste0(feature, ": ", round(value, 1), " nt"),
      TRUE ~ paste0(feature, ": ", round(value, 1))
    )
  }
  
  map_dfr(clusters, function(cl) {
    ids <- df %>% filter(cluster == cl) %>% pull(row_id)
    centroid_x <- median(df$umap_x[df$cluster == cl])
    centroid_y <- median(df$umap_y[df$cluster == cl])
    
    this_mod <- mod_matrix %>%
      filter(row_id %in% ids) %>%
      select(all_of(mod_pos_cols)) %>%
      summarise(across(everything(), mean)) %>%
      pivot_longer(everything(), names_to = "feature", values_to = "this_cluster_value")
    
    this_cont <- df %>%
      filter(cluster == cl) %>%
      summarise(across(all_of(continuous_cols), ~ mean(.x, na.rm = TRUE))) %>%
      pivot_longer(everything(), names_to = "feature", values_to = "this_cluster_value")
    
    this_all <- bind_rows(this_mod, this_cont) %>%
      left_join(overall_all, by = "feature") %>%
      mutate(
        pct_diff = abs(this_cluster_value - overall_value) / pmax(abs(overall_value), 1e-9),
        is_mod = feature %in% mod_pos_cols
      )
    
    # Always include the continuous features, in fixed order
    cont_features <- this_all %>%
      filter(feature %in% continuous_cols) %>%
      arrange(match(feature, continuous_cols))
    
    # Always include the top N mod features by pct_diff
    top_mod_features <- this_all %>%
      filter(is_mod) %>%
      arrange(desc(pct_diff)) %>%
      slice_head(n = top_n_mods)
    
    feature_lines <- bind_rows(cont_features, top_mod_features) %>%
      mutate(line = format_line(feature, this_cluster_value, is_mod)) %>%
      pull(line)
    
    label_text <- paste0(
      "Cluster ", cl, " (n=", length(ids), ")\n",
      paste(feature_lines, collapse = "\n")
    )
    
    tibble(cluster = cl, umap_x = centroid_x, umap_y = centroid_y, label = label_text)
  })
}

cluster_labels <- build_cluster_labels(
  df, mod_matrix,
  continuous_cols = c("read_end", "read_length", "poly_a_length", "read_start"),
  mod_pos_cols,
  top_n_mods = 2
)
# -------------------------------------------------------------------
# 10. Plot: same colors, same layout, same clusters, faceted by file
# -------------------------------------------------------------------
all_cluster_levels <- levels(df$cluster)

df %>%
  filter(cluster != "0") %>%
  ggplot(aes(x = umap_x, y = umap_y, color = cluster)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_discrete(drop = FALSE) +
  facet_wrap(~source_file) +
  labs(title = "UMAP by cluster, split by sample", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

# optional: single combined (non-faceted) plot with labels, for the
# overall cluster structure across all files at once
df %>%
  filter(cluster != "0") %>%
  ggplot(aes(x = umap_x, y = umap_y, color = cluster)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_discrete(drop = FALSE) +
  geom_label_repel(
    data = cluster_labels,
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
  theme_minimal() +
  theme(legend.position = "none")
