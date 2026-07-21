library(tidyverse)
library(umap)
library(randomForest)
library(dbscan)
library(ggrepel)
library(proxy)
library(processx)
# --------------------- FILTER BAMS AND TAKE SCREENSHOTS FOR EACH SAMPLE --------------------- 

bam_files <- c(
  "/Users/jlevendis/01_m6A_3p_readthrough_analysis/01_BAM_filtering_output/28C1_to_pfal.50MAPQ.sorted.bam",
  "/Users/jlevendis/01_m6A_3p_readthrough_analysis/01_BAM_filtering_output/28K1_to_pfal.50MAPQ.sorted.bam"
)

base_dir <- "~/rqc"
igv_path <- "/Applications/IGV_2.19.7.app/Contents/MacOS/IGV"
igv_port <- 60151
genome <- "/Users/jlevendis/Downloads/Pfalciparum3D7/Pfalciparum3D7.genome"
annotation <- "/Users/jlevendis/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"

for(bam_file in bam_files){
  
  sample <- tools::file_path_sans_ext(basename(bam_file))
  message("Processing ", sample)
  
  output_dir <- file.path(base_dir,paste0(sample,"_cluster_bams"))
  snapshot_dir <- file.path(base_dir,paste0(sample,"_igv_screenshots"))
  if(dir.exists(output_dir)) unlink(output_dir, recursive=TRUE)
  if(dir.exists(snapshot_dir)) unlink(snapshot_dir, recursive=TRUE)
  
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(snapshot_dir, recursive=TRUE, showWarnings=FALSE)
  
  clusters <- df %>% filter(cluster!="0") %>% pull(cluster) %>% unique()
  
  for(cl in clusters){
    message("Processing ", cl)
    ids <- df %>% filter(cluster==cl) %>% pull(read_id) %>% unique()
    id_file <- tempfile(fileext=".txt")
    writeLines(ids,id_file)
    
    out_bam <- file.path(output_dir,paste0("cluster_",cl,".bam"))
    
    system(sprintf("samtools view -b -N %s %s > %s",id_file,bam_file,out_bam))
    system(sprintf("samtools index %s",out_bam))
    
    unlink(id_file)
  }
  
  region <- paste0(
    unique(df$contig),":",
    floor(quantile(df$read_start,0.01,na.rm=TRUE)),
    "-",
    ceiling(quantile(df$read_end,0.99,na.rm=TRUE))
  )
  
  igv_process <- process$new(
    igv_path,
    args=c("-p",as.character(igv_port)),
    stdout="|",
    stderr="|"
  )
  
  Sys.sleep(10)
  
  igv_send <- function(cmd){
    con <- socketConnection("localhost",igv_port,blocking=TRUE,open="w")
    writeLines(cmd,con)
    close(con)
    Sys.sleep(1)
  }
  
  igv_send("new")
  igv_send(paste("genome",genome))
  Sys.sleep(5)
  igv_send(paste("load",annotation))
  Sys.sleep(5)
  igv_send(paste("goto",region))
  
  cluster_bams <- list.files(output_dir,pattern="\\.bam$",full.names=TRUE)
  
  cluster_bams <- cluster_bams[order(
    as.numeric(gsub("\\D", "", basename(cluster_bams)))
  )]
  
  for(bam in cluster_bams){
    
    cluster <- tools::file_path_sans_ext(basename(bam))
    message("Loading ", cluster)
    
    igv_send(paste("load", bam))
    Sys.sleep(3)
    
  }
  
  igv_send("squish")
  igv_send("preference SAM.COLOR_BY BASE_MODIFICATION")
  igv_send("preference BASEMOD.THRESHOLD 0.95")
  igv_send("preference SAM.HIDE_SMALL_INDEL TRUE")
  igv_send("preference SAM.HIDE_SMALL_INDEL_BP_THRESHOLD 10")
  Sys.sleep(3)
  
  igv_send(paste("snapshotDirectory",snapshot_dir))
  igv_send(paste("snapshot",paste0(cluster,".png")))
  
  igv_send("exit")
  Sys.sleep(5)
  igv_process$kill()
}