# get in the embryo functions

source("/nfs/research1/marioni/jonny/embryos/scripts/core_functions.R")

load_data = function(
  normalise = TRUE, 
  remove_doublets = FALSE, 
  remove_stripped = FALSE, 
  load_corrected = FALSE,
  remove_pool2 = FALSE){
  
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets = TRUE
    remove_stripped = TRUE
  }
  
  require(scran)
  require(scater)
  require(SingleCellExperiment)
  require(Matrix)
  require(scales)
  
  # counts = readMM("/nfs/research1/marioni/jonny/chimera-t/data/raw_counts.mtx")
  counts = readRDS("/nfs/research1/marioni/jonny/chimera-t/data/raw_counts.rds")
  genes = read.table("/nfs/research1/marioni/jonny/chimera-t/data/genes.tsv", stringsAsFactors = F)
  meta = read.table("/nfs/research1/marioni/jonny/chimera-t/data/meta.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
  
  rownames(counts) = genes[,1] #ensembl
  colnames(counts) = meta$cell
  
  sce = SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs = read.table("/nfs/research1/marioni/jonny/chimera-t/data/sizefactors.tab", stringsAsFactors = F)[,1]
    sizeFactors(sce) = sfs
    sce = scater::normalize(sce)
  }

  if(remove_pool2){
    sce = scater::normalize(sce[, meta$pool != 2])
    meta = meta[meta$pool != 2, ]
  }
  
  if(remove_doublets){
    sce = scater::normalize(sce[,!meta$celltype.mapped == "Doublet"])
    meta = meta[!meta$celltype.mapped == "Doublet",]
  }
  
  if(remove_stripped){
    sce = scater::normalize(sce[,!meta$celltype.mapped == "Stripped"])
    meta = meta[!meta$celltype.mapped == "Stripped", ]
  }
  
  if(load_corrected){
    corrected = readRDS("/nfs/research1/marioni/jonny/chimera-t/data/corrected_pcas.rds")
    assign("corrected", corrected, envir = .GlobalEnv)
    
  }
  
  assign("genes", genes, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  assign("sce", sce, envir = .GlobalEnv)
  
  #add sample colours when meta is available
  small_meta = meta[!duplicated(meta$sample),]
  sample_colours = c(brewer_pal(palette = "Reds")(sum(small_meta$tomato)+1)[-1],
                     brewer_pal(palette = "Greys")(sum(!small_meta$tomato)+1)[-1])
  names(sample_colours) = c(small_meta$sample[small_meta$tomato],
                            small_meta$sample[!small_meta$tomato])
  
  assign("sample_colours", sample_colours, envir = .GlobalEnv)
  
  
  invisible(0)
}


