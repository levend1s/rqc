library(tidyverse)
library(processx)

# --------------------- CONFIG ---------------------
# igv_path   <- "/Applications/IGV_2.19.6.app/Contents/MacOS/IGV"
# genome     <- "/Users/joshualevendis/Documents/RNA/honours/Pfalciparum3D7/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"
# annotation <- "/Users/joshualevendis/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"

labels <- c(
  "28C1", "28K1"
)

# labels <- c(
#   "36C1", "36K1"
# )

igv_path   <- "/Applications/IGV_2.19.8.app/Contents/MacOS/IGV"
genome <- "/Users/jlevendis/Downloads/Pfalciparum3D7/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"
annotation <- "~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"

infile <- c("~/rqc/cluster_transcripts_results.tsv.umap")
df <- read.delim(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Detect which modification columns are present in `df` and store them in
# `mods`. This intersects the expected set with the actual column names.
{
  mod_candidates <- c("m6A", "pseU", "m6A_inosine", "m5C")
  mods <- intersect(mod_candidates, colnames(df))
  message("Detected modification columns: ", if (length(mods)) paste(mods, collapse = ", ") else "(none)")
}

# Define preferred colors for known modifications; absent mods get a grey
# default of 203,203,203
{
  mod_colors <- list(
    m6A = "51,0,111",
    pseU = "253,210,50",
    m6A_inosine = "21,115,17",
    m5C = "255,0,0"
  )

  mod_igv_names <- list(
    m6A = "BASEMOD.A_COLOR",
    pseU = "BASEMOD.OTHER_COLOR",
    m6A_inosine = "BASEMOD.17596_COLOR",
    m5C = "BASEMOD.M_COLOR"
  )
  default_mod_rgb <- "203,203,203"
}

base_dir   <- path.expand("~/rqc/test")

bam_files <- unique(df$bamfile_path[df$label %in% labels])

# IGV screenshot settings
igv_port   <- 60151
MOD_PROB_THRESHOLD <- unique(df$mod_prob_threshold)
stopifnot(length(MOD_PROB_THRESHOLD) == 1)
MOD_PROB_THRESHOLD <- MOD_PROB_THRESHOLD[[1]]
INDEL_THRESHOLD <- 10



# OPTIONS
SKIP_BAM_REGENERATION <- FALSE
DEBUG <- FALSE
# Select the first 17 unique cluster ids from the data frame. If fewer than
# 17 clusters exist, select them all. This replaces the previous "all"
# default so the script processes a manageable subset by default.
all_clusters <- df %>% pull(cluster) %>% as.character() %>% unique()

# CLUSTERS_TO_PROCESS <- as.character(head(all_clusters, 16))
CLUSTERS_TO_PROCESS <- "all"
# CLUSTERS_TO_PROCESS <- setdiff(all_clusters, "cluster2")

# TODO: if many clusters, just process the biggest clusters and print a message


# --------------------- IGV socket helpers ---------------------
# Persistent connection + reading the response is what actually confirms
# IGV finished a command, rather than guessing with Sys.sleep().
igv_connect <- function(port, timeout_s = 30) {
  deadline <- Sys.time() + timeout_s
  repeat {
    con <- tryCatch(
      socketConnection("localhost", port, blocking = TRUE, open = "r+", timeout = 5),
      error = function(e) NULL
    )
    if (!is.null(con)) return(con)
    if (Sys.time() > deadline) stop("Could not connect to IGV on port ", port)
    Sys.sleep(0.5)
  }
}

igv_send <- function(con, cmd, timeout_s = 60) {
  writeLines(cmd, con)
  deadline <- Sys.time() + timeout_s
  resp <- NULL
  while (Sys.time() < deadline) {
    resp <- readLines(con, n = 1)
    if (length(resp) > 0) break
    Sys.sleep(0.1)
  }

  if (grepl("^error", resp, ignore.case = TRUE)) {
    warning("IGV returned an error for '", cmd, "': ", resp)
  }
  return(resp)
}

run_checked <- function(cmd_fmt, ...) {
  cmd <- sprintf(cmd_fmt, ...)
  status <- system(cmd)
  if (status != 0) stop("Command failed (status ", status, "): ", cmd)
}

# --------------------- MAIN LOOP ---------------------
if (!SKIP_BAM_REGENERATION) {
  for (bam_file in bam_files) {
    
    sample <- tools::file_path_sans_ext(basename(bam_file))
    message("Processing ", sample)
    
    output_dir   <- file.path(base_dir, paste0(sample, "_cluster_bams"))
    snapshot_dir <- file.path(base_dir, paste0(sample, "_igv_screenshots"))
    
    if (dir.exists(output_dir))   unlink(output_dir, recursive = TRUE)
    if (dir.exists(snapshot_dir)) unlink(snapshot_dir, recursive = TRUE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)
    
    clusters <- df %>% pull(cluster) %>% as.character() %>% unique()
    
    if (!identical(CLUSTERS_TO_PROCESS, "all")) {
      clusters <- intersect(clusters, as.character(CLUSTERS_TO_PROCESS))
    }
    
    message("Processing clusters: ", paste(clusters, collapse = ", "))
    
    for (cl in clusters) {
      ids <- df %>% filter(cluster == cl & bamfile_path == bam_file) %>% pull(read_id) %>% unique()
      
      # Extract just the "cluster<N>" prefix from a name like "cluster1_blahblah"
      cl_short <- sub("^(cluster[0-9]+).*", "\\1", cl)
      out_bam <- file.path(output_dir, paste0(cl_short, ".bam"))
      
      if (length(ids) == 0) {
        message("Writing empty BAM for cluster ", cl)
        # Header-only BAM (no reads) so downstream tools still get a valid file.
        run_checked("samtools view -b -H %s -o %s", bam_file, out_bam)
        run_checked("samtools index %s", out_bam)
        next
      }
      message("Writing BAM for cluster ", cl, " (", length(ids), " reads)")
      
      id_file <- tempfile(fileext = ".txt")
      writeLines(ids, id_file)
      
      # Extract reads matching ids and sort the resulting BAM before indexing.
      tmp_prefix <- tempfile(pattern = "samtools_sort_")
      run_checked("samtools view -b -N %s %s | samtools sort -T %s -o %s", id_file, bam_file, tmp_prefix, out_bam)
      run_checked("samtools index %s", out_bam)
      unlink(id_file)
    }
  }
}
  






# ------------------- IGV SCREENSHOTS
for (bam_file in bam_files) {
  sample <- tools::file_path_sans_ext(basename(bam_file))
  message("Processing ", sample)
  
  output_dir   <- file.path(base_dir, paste0(sample, "_cluster_bams"))
  snapshot_dir <- file.path(base_dir, paste0(sample, "_igv_screenshots"))
  cluster_bams <- list.files(output_dir, pattern="\\.bam$", full.names=TRUE)
  
  clusters <- df %>% pull(cluster) %>% as.character() %>% unique()
  
  if (!identical(CLUSTERS_TO_PROCESS, "all")) {
    clusters <- intersect(clusters, as.character(CLUSTERS_TO_PROCESS))
  }
  
  message("Processing clusters: ", paste(clusters, collapse = ", "))
  
  # Extract just the "cluster<N>" prefix to match against filenames
  clusters_short <- sub("^(cluster[0-9]+).*", "\\1", clusters)
  
  cluster_bams <- cluster_bams[
    basename(cluster_bams) %in% paste0(clusters_short, ".bam")
  ]
  
  cluster_bams <- cluster_bams[
    order(as.numeric(gsub("\\D","",basename(cluster_bams))))
  ]
  
  if (length(cluster_bams) == 0) {
    warning("No cluster BAMs produced for ", sample, " - skipping IGV step")
    next
  }
  
  if (length(unique(df$contig)) > 1) {
    warning("Multiple contigs in df - using only the first for the region")
  }
  region <- paste0(
    df$contig[1], ":",
    floor(quantile(df$read_start, 0.01, na.rm = TRUE)), "-",
    ceiling(quantile(df$read_end, 0.99, na.rm = TRUE))
  )
  
  # If there are many cluster bams, process them in batches so screenshots
  # contain at most 10 clusters each. Start a fresh IGV session per batch to
  # ensure a clean canvas.
  batch_size <- 10
  if (length(cluster_bams) == 0) next
  batches <- split(cluster_bams, ceiling(seq_along(cluster_bams) / batch_size))

  for (bi in seq_along(batches)) {
    batch_bams <- batches[[bi]]
    message(sprintf("Starting IGV for sample %s (batch %d/%d, %d bams)", sample, bi, length(batches), length(batch_bams)))

    igv_process <- process$new(igv_path, args = c("-p", as.character(igv_port)), stdout = "|", stderr = "|")
    # ensure IGV always gets killed for this batch
    on.exit({
      if (exists("igv_process") && igv_process$is_alive()) igv_process$kill()
    }, add = TRUE)

    Sys.sleep(2)
    if (!igv_process$is_alive()) {
      warning("IGV failed to start for ", sample, " batch ", bi)
      next
    }

    con <- tryCatch(igv_connect(igv_port), error = function(e) NULL)
    if (is.null(con)) {
      warning("Could not connect to IGV for sample ", sample, " batch ", bi)
      if (igv_process$is_alive()) igv_process$kill()
      next
    }

    tryCatch({
      igv_send(con, "new")
      igv_send(con, paste("genome", genome))
      igv_send(con, paste("load", annotation))
      igv_send(con, paste("goto", region))

      for (bam in batch_bams) {
        igv_send(con, paste("load", bam))
      }

      message("making pretty...")
      igv_send(con, "squish")
      igv_send(con, paste("preference BASEMOD.THRESHOLD", MOD_PROB_THRESHOLD))

      for (mod in mod_candidates) {
        if (!mod %in% colnames(df)) {
          attr <- paste("preference", mod_igv_names[[mod]], default_mod_rgb)
          igv_send(con, attr)
        }
        else {
          attr <- paste("preference", mod_igv_names[[mod]], mod_colors[[mod]])
          igv_send(con, attr)
        }
      }
      
      igv_send(con, "preference SAM.COLOR_BY BASE_MODIFICATION")
      igv_send(con, paste("preference SAM.HIDE_SMALL_INDEL_BP_THRESHOLD", INDEL_THRESHOLD))
      igv_send(con, "preference SAM.HIDE_SMALL_INDEL TRUE")
      igv_send(con, "preference SAM.SHOW_SOFT_CLIPPED TRUE")

      if (DEBUG) {
        readline(sprintf("Batch %d/%d ready. Press [Enter] to continue (Ctrl+C to interrupt): ", bi, length(batches)))
      }

      igv_send(con, paste("snapshotDirectory", snapshot_dir))
      igv_send(con, paste("snapshot", paste0(sample, "_clusters_batch_", bi, ".png")))
    }, error = function(e) {
      warning("IGV batch ", bi, " failed: ", conditionMessage(e))
    }, finally = {
      try(close(con), silent = TRUE)
      if (igv_process$is_alive()) igv_process$kill()
      Sys.sleep(1)
    })
  }
}