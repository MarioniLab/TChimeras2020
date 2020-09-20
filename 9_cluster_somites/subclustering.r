library(scran)
library(WGCNA)
library(igraph)
library(biomaRt)
library(pheatmap)

source("/nfs/research1/marioni/jonny/embryos/scripts/core_functions.R")

# Read Atlas
load_data(normalise = TRUE,
  remove_doublets = TRUE,
  remove_stripped = TRUE, 
  load_corrected  = FALSE)

##############
############## 
# Subsetting, HVGs and PCA

ct1 <- gsub(" ",".",meta$celltype) == "Paraxial.mesoderm" & meta$stage == "E8.5"
ct2 <- gsub(" ",".",meta$celltype) == "Somitic.mesoderm" & meta$stage == "E8.5"

sce_pm  <- sce[,ct1]
sce_sm  <- sce[,ct2]

hvgs_pm <- getHVGs(sce_pm)
hvgs_sm <- getHVGs(sce_sm)

pca_pm <- prcomp(t(logcounts(sce_pm)[hvgs_pm,]))
pca_sm <- prcomp(t(logcounts(sce_sm)[hvgs_sm,]))

##############
##############
# Correct pc's for batch effect

index  <- match(colnames(sce_pm), meta$cell)

order_df = meta[index,][!duplicated(meta$sample[index]), c("stage", "sample")]
order_df$ncells = sapply(order_df$sample, function(x) sum(meta$sample[index] == x))
order_df$stage  = factor(order_df$stage,  levels = rev(c("E8.5")))
  order_df = order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage = as.character(order_df$stage)

pca_corrected_pm <- doBatchCorrect(counts = logcounts(sce_pm)[hvgs_pm, ], 
    timepoints      = meta$stage[index], 
    samples         = meta$sample[index], 
    timepoint_order = order_df$stage, 
    sample_order    = order_df$sample, 
    pc_override     = pca_pm$x[,1:50])


index  <- match(colnames(sce_sm), meta$cell)

order_df = meta[index,][!duplicated(meta$sample[index]), c("stage", "sample")]
order_df$ncells = sapply(order_df$sample, function(x) sum(meta$sample[index] == x))
order_df$stage  = factor(order_df$stage, levels = rev(c("E8.5")))
  order_df = order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage = as.character(order_df$stage)

pca_corrected_sm <- doBatchCorrect(counts = logcounts(sce_sm)[hvgs_sm, ], 
    timepoints      = meta$stage[index], 
    samples         = meta$sample[index], 
    timepoint_order = order_df$stage, 
    sample_order    = order_df$sample, 
    pc_override     = pca_sm$x[,1:50])

##############
##############
# Louvain clustering / annotation

graph_pm <- buildSNNGraph(t(pca_corrected_pm), k=15, d=NA)

graph_sm <- buildSNNGraph(t(pca_corrected_sm), k=15, d=NA)

clusters_pm <- igraph::cluster_louvain(graph_pm)$membership
names(clusters_pm) <- colnames(sce_pm)
table(clusters_pm)

clusters_sm <- igraph::cluster_louvain(graph_sm)$membership
names(clusters_sm) <- colnames(sce_sm)
table(clusters_sm)

colours_pm <- labels2colors(clusters_pm)
colours_sm <- labels2colors(clusters_sm)
names(colours_pm) <- names(clusters_pm)
names(colours_sm) <- names(clusters_sm)

subtypes_pm <- colours_pm
subtypes_pm <- gsub("brown",     "Sclerotome", subtypes_pm)
subtypes_pm <- gsub("blue",      "Anterior-most_somites", subtypes_pm)
subtypes_pm <- gsub("yellow",    "Dermomyotome", subtypes_pm)
subtypes_pm <- gsub("turquoise", "Head_mesoderm", subtypes_pm)
colours_annot_pm <- labels2colors(subtypes_pm, 
                    colorSeq=c("red", "brown", "darkgreen", "yellow")) 

subtypes_sm <- colours_sm
subtypes_sm <- gsub("brown",     "Posterior-most_somites", subtypes_sm)
subtypes_sm <- gsub("turquoise", "Presomitic_mesoderm", subtypes_sm)
subtypes_sm <- gsub("blue",      "Presomitic_mesoderm", subtypes_sm)
colours_annot_sm <- labels2colors(subtypes_sm)

paraxialmes_meta <- data.frame(
  cell = names(colours_pm),
  louvain_cluster = clusters_pm,
  louvain_colours = colours_pm, 
  celltype = subtypes_pm, 
  celltype_colours = colours_annot_pm) 

somiticmes_meta  <- data.frame(
  cell = names(colours_sm),
  louvain_cluster = clusters_sm,
  louvain_colours = colours_sm, 
  celltype = subtypes_sm, 
  celltype_colours = colours_annot_sm) 

write.table(paraxialmes_meta,
 file = "paraxialmeso_e85_annotation.txt", row.names=FALSE, quote=FALSE)

write.table(somiticmes_meta,
 file = "somiticmeso_e85_annotation.txt", row.names=FALSE, quote=FALSE)

##############
##############
# Cell sets (WOT)
 
cellSet <- NULL
for (i in 1:length(unique(subtypes_pm))){
 celltype_index <- which(unique(subtypes_pm)[i] == subtypes_pm)
 cell_set     <- rep(NA,max(table(subtypes_pm))+2)
 cell_set_tmp <- c(unique(subtypes_pm)[i],"cell_type", names(colours_pm)[celltype_index])
 cell_set[1:length(cell_set_tmp)]  <- cell_set_tmp
 cell_set <- gsub(" ",".",cell_set)
 cellSet   <- rbind(cellSet, cell_set)  
}

write.table(cellSet, file=paste0(path2data,"cell_sets_paraxialmes_subclustering_wot.gmt"),
  row.names=FALSE, col.names=FALSE, sep="\t", na="", quote=FALSE)

cellSet <- NULL
for (i in 1:length(unique(subtypes_sm))){
 celltype_index <- which(unique(subtypes_sm)[i] == subtypes_sm)
 cell_set     <- rep(NA,max(table(subtypes_sm))+2)
 cell_set_tmp <- c(unique(subtypes_sm)[i],"cell_type", names(colours_sm)[celltype_index])
 cell_set[1:length(cell_set_tmp)]  <- cell_set_tmp
 cell_set <- gsub(" ",".",cell_set)
 cellSet   <- rbind(cellSet, cell_set)  
}

write.table(cellSet, file=paste0(path2data,"cell_sets_somiticmes_subclustering_wot.gmt"),
  row.names=FALSE, col.names=FALSE, sep="\t", na="", quote=FALSE)

##############
##############
# Find markers - heatmap
 
paraxialmes_meta <- read.table(file = "paraxialmeso_e85_annotation.txt", header = TRUE)

somiticmes_meta  <- read.table(file = "somiticmeso_e85_annotation.txt", header = TRUE)

meta_chimera <- rbind(somiticmes_meta, paraxialmes_meta)

clus_markers <- findMarkers(cbind(
  sce_atlas_sm[union(hvgs_sm, hvgs_pm),],
  sce_atlas_pm[union(hvgs_sm, hvgs_pm),]), 
  meta_chimera$celltype, assay.type = "logcounts", get.spikes = FALSE)
  
cluster_names <- c("Presomitic_mesoderm", "Posterior-most_somites", 
  "Head_mesoderm", "Dermomyotome", "Sclerotome", "Anterior-most_somites")

genes <- NULL
for (cl in cluster_names){
  markers     <- clus_markers[[cl]]
  top.markers <- rownames(marker)[marker$Top <= 10]
  genes <- unique(c(genes, top.markers))
}
 
genes <- c(genes,
  "ENSMUSG00000000567",
  "ENSMUSG00000001497",
  "ENSMUSG00000046714",
  "ENSMUSG00000004872",
  "ENSMUSG00000000435",
  "ENSMUSG00000051367",
  "ENSMUSG00000028023")

data.mat <- cbind(logcounts(sce_atlas_sm)[genes,], logcounts(sce_atlas_pm)[genes,])
clus.ann <- meta_chimera$celltype

means.table <- cbind(
  cl1 = apply(data.mat[, which(clus.ann == cluster_names[1])], 1, mean),
  cl2 = apply(data.mat[, which(clus.ann == cluster_names[2])], 1, mean),
  cl3 = apply(data.mat[, which(clus.ann == cluster_names[3])], 1, mean),
  cl4 = apply(data.mat[, which(clus.ann == cluster_names[4])], 1, mean),
  cl5 = apply(data.mat[, which(clus.ann == cluster_names[5])], 1, mean),
  cl6 = apply(data.mat[, which(clus.ann == cluster_names[6])], 1, mean))
colnames(means.table) <- cluster_names

mouse_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_map <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),
            filters = "ensembl_gene_id", values = rownames(means.table), 
            mart = mouse_ensembl)
genes    <- gene_map$mgi_symbol[match(rownames(means.table), gene_map$ensembl_gene_id)]

rownames(means.table) <- genes

means.table <- means.table[,
  c("Head_mesoderm",
  "Anterior-most_somites",
  "Sclerotome", "Dermomyotome",
  "Posterior-most_somites",
  "Presomitic_mesoderm")]
  
pheatmap(means.table, scale = "row", cellwidth = 15, cellheight = 12, fontsize = 8,
  cluster_cols = FALSE)

