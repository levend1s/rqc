library(tidyverse)
library(processx)

# --------------------- CONFIG ---------------------
# bam_files <- c(
#   "/Users/joshualevendis/Downloads/bams/28C1_to_pfal.50MAPQ.sorted.bam",
#   "/Users/joshualevendis/Downloads/bams/28K1_to_pfal.50MAPQ.sorted.bam"
# )

# igv_path   <- "/Applications/IGV_2.19.6.app/Contents/MacOS/IGV"
# genome     <- "/Users/joshualevendis/Documents/RNA/honours/Pfalciparum3D7/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"
# annotation <- "/Users/joshualevendis/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"

bam_files <- c(
  "/Users/jlevendis/01_m6A_3p_readthrough_analysis/01_BAM_filtering_output/28C1_to_pfal.50MAPQ.sorted.bam",
  "/Users/jlevendis/01_m6A_3p_readthrough_analysis/01_BAM_filtering_output/28K1_to_pfal.50MAPQ.sorted.bam"
)
# bam_files <- unique(df$bamfile_path)

igv_path   <- "/Applications/IGV_2.19.7.app/Contents/MacOS/IGV"
genome <- "/Users/jlevendis/Downloads/Pfalciparum3D7/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"
annotation <- "~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"
igv_port   <- 60151
base_dir   <- path.expand("~/rqc")
base_dir   <- path.expand("~/rqc/test")

CLUSTERS_TO_PROCESS <- "all"  # or e.g. c("1","3","5")
# CLUSTERS_TO_PROCESS <- c("13")

stopifnot(exists("df"))  # this script depends on `df` from the clustering script -
# fail fast and loudly rather than silently using a stale df

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
  if (length(resp) == 0) stop("IGV did not respond to: ", cmd)
  if (grepl("^error", resp, ignore.case = TRUE)) {
    warning("IGV returned an error for '", cmd, "': ", resp)
  }
  invisible(resp)
}

run_checked <- function(cmd_fmt, ...) {
  cmd <- sprintf(cmd_fmt, ...)
  status <- system(cmd)
  if (status != 0) stop("Command failed (status ", status, "): ", cmd)
}

# --------------------- MAIN LOOP ---------------------
for (bam_file in bam_files) {
  
  sample <- tools::file_path_sans_ext(basename(bam_file))
  message("Processing ", sample)
  
  output_dir   <- file.path(base_dir, paste0(sample, "_cluster_bams"))
  snapshot_dir <- file.path(base_dir, paste0(sample, "_igv_screenshots"))
  if (dir.exists(output_dir))   unlink(output_dir, recursive = TRUE)
  if (dir.exists(snapshot_dir)) unlink(snapshot_dir, recursive = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # clusters <- df %>% filter(cluster != "0") %>% pull(cluster) %>% unique()
  clusters <- df %>% pull(cluster) %>% as.character() %>% unique()
  
  if(!identical(CLUSTERS_TO_PROCESS,"all")){
    clusters <- intersect(clusters, as.character(CLUSTERS_TO_PROCESS))
  }
  
  message("Processing clusters: ", paste(clusters, collapse=", "))
  
  for (cl in clusters) {
    ids <- df %>% filter(cluster == cl) %>% pull(read_id) %>% unique()
    if (length(ids) == 0) {
      message("Skipping empty cluster ", cl)
      next
    }
    message("Writing BAM for cluster ", cl, " (", length(ids), " reads)")
    
    id_file <- tempfile(fileext = ".txt")
    writeLines(ids, id_file)
    out_bam <- file.path(output_dir, paste0("cluster_", cl, ".bam"))
    
    run_checked("samtools view -b -N %s %s > %s", id_file, bam_file, out_bam)
    run_checked("samtools index %s", out_bam)
    unlink(id_file)
  }
  
  cluster_bams <- list.files(output_dir, pattern="\\.bam$", full.names=TRUE)
  
  cluster_bams <- cluster_bams[
    basename(cluster_bams) %in% paste0("cluster_",clusters,".bam")
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
  
  igv_process <- process$new(igv_path, args = c("-p", as.character(igv_port)),
                             stdout = "|", stderr = "|")
  
  # ensure IGV always gets killed, even if something below errors
  on.exit({
    if (igv_process$is_alive()) igv_process$kill()
  }, add = TRUE)
  
  # wait for the IGV process itself to be alive before even trying the socket
  Sys.sleep(2)
  if (!igv_process$is_alive()) stop("IGV failed to start for ", sample)
  
  con <- igv_connect(igv_port)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  
  igv_send(con, "new")
  igv_send(con, paste("genome", genome))
  igv_send(con, paste("load", annotation))
  igv_send(con, paste("goto", region))
  
  for (bam in cluster_bams) {
    message("Loading ", basename(bam))
    igv_send(con, paste("load", bam))
  }
  
  message("making pretty...")
  igv_send(con, "squish")
  igv_send(con, "preference SAM.COLOR_BY BASE_MODIFICATION")
  igv_send(con, "preference BASEMOD.THRESHOLD 0.95")
  igv_send(con, "preference SAM.HIDE_SMALL_INDEL TRUE")
  igv_send(con, "preference SAM.HIDE_SMALL_INDEL_BP_THRESHOLD 10")
  # igv_send(con, "preference SAM.SHOW_SOFT_CLIPPED TRUE")
  
  # "OK" from these commands confirms IGV finished re-rendering with the
  # new display settings, so a screenshot here reflects the final state -
  # no fixed-length guess-sleep required
  igv_send(con, paste("snapshotDirectory", snapshot_dir))
  igv_send(con, paste("snapshot", paste0(sample, "_all_clusters.png")))
  
  # igv_send(con, "exit")
  try(close(con), silent = TRUE)
  igv_process$kill()
  Sys.sleep(1)
}