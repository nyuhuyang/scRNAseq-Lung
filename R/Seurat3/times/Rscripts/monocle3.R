####################################
library(monocle3)
library(Seurat)
library(dplyr)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
set.seed=101
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# samples
samples = c("combined","distal","proximal","terminal")
(con <- samples[args])


# 5.0 Preliminaries: Load the data
(load(file = paste0("data/Lung_23",con,"_20190824.Rda")))
Idents(object) = "cell.type"
Epi <- subset(object, idents=c("Alveolar cells/Distal secretory cells",
                               "Ciliated cells","Basal cells","Secretory cells"))
Epi@meta.data = cbind(Epi@meta.data,Epi@reductions$umap@cell.embeddings)
UMAPPlot.1(Epi,cols = ExtractMetaColor(Epi),label.repel = T,label = T,no.legend = T,
           do.print = T, title = "UMAP plot for Epithelial cells",unique.name = "conditions")
if(con =="combined"){
    Epi <- subset(x = Epi, subset = UMAP_1 > -5)
    Epi <- subset(x = Epi, subset = UMAP_2 > 0)
}

if(con =="proximal"){
    Epi <- subset(x = Epi, subset = UMAP_1 <5)
    Epi <- subset(x = Epi, subset = UMAP_2 <0)
}
if(con =="distal"){
    Epi <- subset(x = Epi, subset = UMAP_1 <0)
    Epi <- subset(x = Epi, subset = UMAP_2 > -5)
}
if(con =="terminal"){
    Epi <- subset(x = Epi, subset = UMAP_1 > -5)
    Epi <- subset(x = Epi, subset = UMAP_2 > 0)
}
object <- Epi
remove(Epi);GC()
#Construct monocle cds
cds <- new_cell_data_set(object@assays$SCT@counts,cell_metadata = object@meta.data)
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100, residual_model_formula_str = "~ orig.ident")
table(rownames(reducedDims(cds)[["PCA"]]) == rownames(object@reductions$pca@cell.embeddings))
reducedDims(cds)[["PCA"]] = object@reductions$pca@cell.embeddings
colnames(reducedDims(cds)[["PCA"]]) %<>% gsub("_","",.)
## Step 2: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)
reducedDims(cds)$UMAP <- object@reductions$umap@cell.embeddings
colnames(reducedDims(cds)$UMAP) = NULL
## Step 3: Cluster the cells
cds <- cluster_cells(cds)

## Step 4: Learn a graph
cds <- learn_graph(cds)

## Step 5: Order cells
#cds <- order_cells(cds)

g <- plot_cells(cds, color_cells_by="cell.type",cell_size = 1,
                show_trajectory_graph = TRUE,trajectory_graph_segment_size = 0.75,
                label_cell_groups = F,
                label_groups_by_cluster = F, 
                label_branch_points = F)+scale_colour_manual(values = ExtractMetaColor(object))
jpeg(paste0(path,con,"_pseudotime.jpeg"), units="in", width=10, height=7,res=600)
print(g)+ ggtitle(paste("Pseudotime trajectory of",con,"samples"))+ 
        theme(plot.title = element_text(size=20, hjust = 0.5,face="plain"))
dev.off()