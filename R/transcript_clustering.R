library(tidyverse)
library(umap)
library(randomForest)
library(dbscan)
library(ggrepel)
library(proxy)
library(processx)
library(edgeR)

all_continuous_cols <- c(
  "read_start",
  "read_end",
  "read_length",
  "poly_a_length"
  # "average_quality"
)

# also used for dropping rare features
NUM_SAMPLES <- 4
MIN_CLUSTER_SIZE <- 10 * NUM_SAMPLES
MAX_CLUSTER_SIZE <- 200 * NUM_SAMPLES
MIN_CLUSTER_PERC <- 0.01

UMAP_MIN_DIST <- 0.3
manual_lib_sizes <- c(2448848, 1350852, 1790844, 2283056)  # example values
USE_RF <- FALSE

MOD_TYPES <- c("m6A", "m5C", "pseU", "m6A_inosine")
MOD_TYPES <- c("m6A")
# MOD_TYPES <- NULL

# balanced
w_mod <- 0.3
if (is.null(MOD_TYPES)) w_mod <- 0 else
w_intron <- 0.3
w_cont <- 0.4

# no cont variables
# w_mod <- 0.5
# w_intron <- 0.5
# w_cont <- 0

set.seed(42)

file <- c("~/rqc/cluster_transcripts_results.tsv")

# -------------------------------------------------------------------
# 8. Cluster summaries and labels
# -------------------------------------------------------------------
build_cluster_labels <- function(df,
                                 mod_matrix,
                                 intron_matrix,
                                 continuous_cols,
                                 mod_pos_cols,
                                 intron_cols,
                                 top_n=5){
  
  feature_cols <- c(mod_pos_cols, intron_cols)
  
  if(USE_RF){
    feature_df <- bind_cols(
      df %>% select(all_of(continuous_cols)),
      mod_matrix %>% select(all_of(mod_pos_cols)),
      intron_matrix %>% select(all_of(intron_cols))
    )
    
    rf <- randomForest(
      x=feature_df,
      y=df$cluster,
      importance=TRUE
    )
    
    importance_df <- tibble(
      feature=rownames(importance(rf)),
      loading=importance(rf)[,"MeanDecreaseGini"]
    )
  } else{
    importance_df <- tibble(
      feature=c(continuous_cols, feature_cols),
      loading=1
    )
  }
  
  overall_cont <- df %>%
    summarise(across(all_of(continuous_cols), ~mean(.x, na.rm=TRUE))) %>%
    pivot_longer(everything(), names_to="feature", values_to="overall")
  
  overall_mod <- if(length(mod_pos_cols)>0){
    mod_matrix %>%
      select(all_of(mod_pos_cols)) %>%
      summarise(across(everything(), ~mean(.x, na.rm=TRUE))) %>%
      pivot_longer(everything(), names_to="feature", values_to="overall")
  } else{
    tibble(feature=character(), overall=numeric())
  }
  
  overall_intron <- if(length(intron_cols)>0){
    intron_matrix %>%
      select(all_of(intron_cols)) %>%
      summarise(across(everything(), ~mean(.x, na.rm=TRUE))) %>%
      pivot_longer(everything(), names_to="feature", values_to="overall")
  } else{
    tibble(feature=character(), overall=numeric())
  }
  
  overall <- bind_rows(overall_cont, overall_mod, overall_intron)
  
  feature_sd <- bind_cols(
    df %>% select(all_of(continuous_cols)),
    mod_matrix %>% select(all_of(mod_pos_cols)),
    intron_matrix %>% select(all_of(intron_cols))
  ) %>%
    summarise(across(everything(), ~sd(.x, na.rm=TRUE))) %>%
    pivot_longer(everything(), names_to="feature", values_to="sd")
  
  clusters <- setdiff(levels(df$cluster), "0")
  
  map_dfr(clusters, function(cl){
    
    ids <- df %>%
      filter(cluster==cl) %>%
      pull(read_id)
    
    cont_summary <- df %>%
      filter(cluster==cl) %>%
      summarise(across(all_of(continuous_cols), ~mean(.x, na.rm=TRUE))) %>%
      pivot_longer(everything(), names_to="feature", values_to="value")
    
    mod_summary <- if(length(mod_pos_cols)>0){
      mod_matrix %>%
        filter(read_id %in% ids) %>%
        select(all_of(mod_pos_cols)) %>%
        summarise(across(everything(), ~mean(.x, na.rm=TRUE))) %>%
        pivot_longer(everything(), names_to="feature", values_to="value")
    } else{
      tibble(feature=character(), value=numeric())
    }
    
    intron_summary <- if(length(intron_cols)>0){
      intron_matrix %>%
        filter(read_id %in% ids) %>%
        select(all_of(intron_cols)) %>%
        summarise(across(everything(), ~mean(.x, na.rm=TRUE))) %>%
        pivot_longer(everything(), names_to="feature", values_to="value")
    } else{
      tibble(feature=character(), value=numeric())
    }
    
    summary <- bind_rows(cont_summary, mod_summary, intron_summary) %>%
      left_join(overall, by="feature") %>%
      left_join(feature_sd, by="feature") %>%
      left_join(importance_df, by="feature") %>%
      mutate(
        sd=replace_na(sd,1),
        loading=replace_na(loading,1),
        effect=abs(value-overall)/pmax(sd,1e-8),
        score=effect*loading
      ) %>%
      arrange(desc(score))
    
    top_features <- head(summary, top_n)
    
    feature_text <- map_chr(seq_len(nrow(top_features)), function(i){
      
      f <- top_features$feature[i]
      v <- top_features$value[i]
      l <- top_features$loading[i]
      
      if(f %in% mod_pos_cols){
        txt <- sprintf("%s %.0f%%",f,100*v)
      } else if(f %in% intron_cols){
        txt <- sprintf("Intron %s %.0f%%",sub("^intron_","",f),100*v)
      } else if(f %in% c("read_start","read_end","read_length")){
        txt <- sprintf("%s %.0f bp",f,v)
      } else if(f=="poly_a_length"){
        txt <- sprintf("%s %.1f nt",f,v)
      } else if(f=="average_quality"){
        txt <- sprintf("Q %.1f",v)
      } else{
        txt <- sprintf("%s %.2f",f,v)
      }
      
      if(USE_RF){
        sprintf("%s (RF %.1f)",txt,l)
      } else{
        txt
      }
    })
    
    tibble(
      cluster=cl,
      n=length(ids),
      umap_x=median(df$umap_x[df$cluster==cl]),
      umap_y=median(df$umap_y[df$cluster==cl]),
      label=paste0(
        "Cluster ",cl,
        "\n(n=",length(ids),")",
        "\n",
        paste(feature_text,collapse="\n")
      )
    )
  })
}

#---------------------------

df <- read.delim(file)

make_mod_matrix <- function(df, mod_type){
  col <- paste0(mod_type, "_positions")
  
  df_long <- df %>%
    select(read_id, all_of(col)) %>%
    mutate(
      !!col := str_remove_all(.data[[col]], "\\[|\\]"),
      !!col := na_if(.data[[col]], ""),
      !!col := str_split(.data[[col]], ",\\s*")
    ) %>%
    unnest(all_of(col)) %>%
    filter(!is.na(.data[[col]]), .data[[col]] != "") %>%
    rename(mod_position = all_of(col)) %>%
    count(read_id, mod_position, name="count")
  
  df %>%
    select(read_id) %>%
    left_join(
      df_long %>%
        pivot_wider(
          names_from=mod_position,
          values_from=count,
          values_fill=0,
          names_prefix=paste0(mod_type,"_")
        ),
      by="read_id"
    ) %>%
    replace(is.na(.),0)
}

mod_matrices <- lapply(
  MOD_TYPES,
  function(x) make_mod_matrix(df, x)
)

names(mod_matrices) <- MOD_TYPES

mod_matrix <- df %>%
  select(read_id)

for(i in seq_along(mod_matrices)){
  mod_matrix <- mod_matrix %>%
    left_join(mod_matrices[[i]], by="read_id")
}

mod_matrix[is.na(mod_matrix)] <- 0

all_mod_pos_cols <- grep(
  "^m6A_|^m5C_|^pseU_|^m6A_inosine_",
  names(mod_matrix),
  value=TRUE
)

mod_X <- mod_matrix %>%
  select(all_of(all_mod_pos_cols)) %>%
  as.matrix()

# make intron matrix
make_intron_matrix <- function(df){
  
  df_long <- df %>%
    select(read_id, introns) %>%
    mutate(
      introns = str_remove_all(introns, "\\[|\\]"),
      introns = na_if(introns, ""),
      introns = str_split(introns, "\\),\\s*\\(")
    ) %>%
    mutate(
      introns = lapply(
        introns,
        function(x){
          gsub("\\(|\\)", "", x)
        }
      )
    ) %>%
    unnest(introns) %>%
    filter(!is.na(introns), introns != "") %>%
    mutate(
      intron_id = paste0("intron_", introns)
    ) %>%
    distinct(read_id, intron_id) %>%
    mutate(value=1)
  
  df %>%
    select(read_id) %>%
    left_join(
      df_long %>%
        select(read_id, intron_id, value) %>%
        pivot_wider(
          names_from=intron_id,
          values_from=value,
          values_fill=0
        ),
      by="read_id"
    ) %>%
    replace(is.na(.),0)
}

intron_matrix <- make_intron_matrix(df)

intron_X <- intron_matrix %>%
  select(-read_id) %>%
  as.matrix()

# DEBUG
# intron_counts <- colSums(intron_matrix[-1] > 0)
# 
# intron_summary <- tibble(
#   intron = names(intron_counts),
#   n_reads = intron_counts,
#   pct_reads = 100 * intron_counts / nrow(intron_matrix)
# ) %>%
#   arrange(desc(n_reads))
# 
# print(intron_summary, n = 50)

# FILTERING
# mod columns
min_cluster_size <- min(MAX_CLUSTER_SIZE, max(MIN_CLUSTER_SIZE, ceiling(MIN_CLUSTER_PERC * nrow(df) / NUM_SAMPLES)))

mod_keep <- colSums(mod_X > 0) >= min_cluster_size
mod_X <- mod_X[, mod_keep, drop=FALSE]

# intron columns
intron_keep <- colSums(intron_X > 0) >= min_cluster_size
intron_X <- intron_X[, intron_keep, drop=FALSE]

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
intron_X_tfidf <- compute_freq_weighted(intron_X)

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
mod_tfidf_cosine_dist <- as.matrix(proxy::dist(mod_X_tfidf, method = "cosine"))
mod_tfidf_cosine_dist[is.na(mod_tfidf_cosine_dist)] <- 1

intron_tfidf_cosine_dist <- as.matrix(proxy::dist(intron_X_tfidf, method = "cosine"))
intron_tfidf_cosine_dist[is.na(intron_tfidf_cosine_dist)] <- 1

euclidean_dist <- as.matrix(dist(cont_X, method = "euclidean"))

normalize01 <- function(m) (m - min(m)) / (max(m) - min(m))

mod_tfidf_cosine_dist_norm <- normalize01(mod_tfidf_cosine_dist)
intron_tfidf_cosine_dist_norm <- normalize01(intron_tfidf_cosine_dist)
euclidean_dist_norm <- normalize01(euclidean_dist)

combined_dist <- w_mod * mod_tfidf_cosine_dist_norm + w_cont * euclidean_dist_norm + w_intron * intron_tfidf_cosine_dist_norm

# -------------------------------------------------------------------
# 6. Run UMAP ONCE on the combined data -> shared coordinate space
# -------------------------------------------------------------------
# TODO: this doesn't normalise by sample variation

custom_config <- umap.defaults
custom_config$input <- "dist"
custom_config$n_neighbors <- min_cluster_size
custom_config$min_dist <- UMAP_MIN_DIST
umap_result <- umap(combined_dist, config = custom_config)

df$umap_x <- umap_result$layout[, 1]
df$umap_y <- umap_result$layout[, 2]

# -------------------------------------------------------------------
# 7. Cluster ONCE on the combined distance matrix -> shared cluster IDs
# -------------------------------------------------------------------
clust <- hdbscan(umap_result$layout, minPts = min_cluster_size)

df$cluster <- as.factor(clust$cluster)
table(df$cluster, df$label)  # check cluster x file distribution -- see batch-effect note below

# -------------------------------------------------------------------
# 9. Plot: same colors, same layout, same clusters, faceted by file
# -------------------------------------------------------------------
cluster_labels <- build_cluster_labels(
  df,
  mod_matrix,
  intron_matrix,
  continuous_cols,
  colnames(mod_X),
  colnames(intron_X),
  top_n = 5
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
# 10. Differential abundance of clusters across samples (edgeR)
# -------------------------------------------------------------------
# TODO stop ignoring cluster 0!
# clusters x samples count matrix
cluster_counts <- table(df$cluster, df$label)
cluster_counts <- as.matrix(cluster_counts)

# drop the noise cluster ("0" from HDBSCAN) before testing
# cluster_counts <- cluster_counts[rownames(cluster_counts) != "0", ]

# metadata: map label -> condition/group
# (edit this to reflect your actual experimental design)
sample_info <- data.frame(
  label = colnames(cluster_counts),
  group = c("control", "control", "treatment", "treatment") # <- your groups
)

y <- DGEList(counts = cluster_counts, group = sample_info$group, lib.size = manual_lib_sizes)

# library size = total reads per sample, TMM normalizes for
# compositional differences (accounts for the fact that if one
# cluster balloons in a sample, it mechanically shrinks others)
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ group, data = sample_info)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)  # tests conditionB vs conditionA

topTags(qlf, n = Inf)

volcano_df <- topTags(qlf, n = Inf)$table %>%
  as_tibble(rownames = "cluster") %>%
  mutate(
    sig = case_when(
      FDR < 0.05 & logFC >  1 ~ "Up",
      FDR < 0.05 & logFC < -1 ~ "Down",
      TRUE                    ~ "NS"
    )
  )

ggplot(volcano_df, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_label_repel(
    data = volcano_df %>% filter(sig != "NS"),
    aes(label = paste0("Cluster ", cluster)),
    size = 3,
    fontface = "bold",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "NS" = "grey70")) +
  labs(
    title = "Differential cluster abundance",
    subtitle = "logFC and significance from edgeR QLF test",
    x = "log2 fold change (treatment / control)",
    y = expression(-log[10](p-value)),
    color = "Significance"
  ) +
  theme_classic()

