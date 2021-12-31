# conda activate r4.1.1
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(data.table)
set.seed(1234)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#=========== sample 30 ==================
object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
D <- subset(object, subset = Family %in% c("ASE","AT")
               &  Regions == "distal")
COPD <- subset(object, subset = Family %in% c("ASE","AT")
                   &  Regions == "COPD")
# Building trajectories with Monocle 3
D.cds <- as.cell_data_set(D)
D.cds <- cluster_cells(cds = D.cds, reduction_method = "UMAP")
D.cds <- learn_graph(D.cds, use_partition = TRUE)
D.cds <- order_cells(D.cds)


COPD.cds <- as.cell_data_set(COPD)
COPD.cds <- cluster_cells(cds = COPD.cds, reduction_method = "UMAP")
COPD.cds <- learn_graph(COPD.cds, use_partition = TRUE)
COPD.cds <- order_cells(COPD.cds)

# plot trajectories colored by pseudotime
jpeg(paste0(path,"distal_Cell_subtype.jpeg"), units="in", width=7, height=7,res=600)
plot_cells(
    cds = D.cds,
    color_cells_by = "Cell_subtype",
    show_trajectory_graph = TRUE
)
dev.off()

jpeg(paste0(path,"distal_pseudotime.jpeg"), units="in", width=7, height=7,res=600)
plot_cells(
    cds = D.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
)
dev.off()

jpeg(paste0(path,"COPD_Cell_subtype.jpeg"), units="in", width=7, height=7,res=600)
plot_cells(
    cds = COPD.cds,
    color_cells_by = "Cell_subtype",
    show_trajectory_graph = TRUE
)
dev.off()

jpeg(paste0(path,"COPD_pseudotime1.jpeg"), units="in", width=7, height=7,res=600)
plot_cells(
    cds = COPD.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
)
dev.off()

D <- AddMetaData(object = D, 
                 metadata = D.cds@principal_graph_aux@listData$UMAP$pseudotime,
                 col.name = "pseudotime")
D@meta.data %<>% cbind(D[["umap"]]@cell.embeddings)
D$barcode = colnames(D)
COPD <- AddMetaData(object = COPD, 
                 metadata = COPD.cds@principal_graph_aux@listData$UMAP$pseudotime,
                 col.name = "pseudotime")
COPD@meta.data %<>% cbind(COPD[["umap"]]@cell.embeddings)
COPD$barcode = colnames(COPD)

fwrite(D@meta.data[,c("barcode","pseudotime","UMAP_1","UMAP_2")],
       file = paste0(path,"distal_pseudotime.csv"))

fwrite(COPD@meta.data[,c("barcode","pseudotime","UMAP_1","UMAP_2")],
       file = paste0(path,"COPD_pseudotime.csv"))
