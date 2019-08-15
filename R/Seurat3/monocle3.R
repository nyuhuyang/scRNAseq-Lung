####################################
# install packages
BiocManager::install()
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment'))
install.packages("reticulate")
reticulate::py_install("louvain")
devtools::install_github('cole-trapnell-lab/monocle3')
####################################
library(monocle3)
library(Seurat)
library(dplyr)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 5.0 Preliminaries: Load the data

(load(file = "data/Lung_9_20190813.Rda"))
Idents(object) = "cell.type"
UMAPPlot.1(object,cols = ExtractMetaColor(object),label.repel = T,label = T,
           no.legend = T,do.print = T)
Epi <- subset(object, idents=c("Alveolar cells/Distal secretory cells",
                               "Ciliated cells","Basal cells","Secretory cells"))
Epi@meta.data = cbind(Epi@meta.data,Epi@reductions$umap@cell.embeddings)
UMAPPlot.1(Epi,cols = ExtractMetaColor(Epi),label.repel = T,label = T,no.legend = T,
           do.print = T, title = "UMAP plot for Epithelial cells")
Epi <- subset(x = Epi, subset = UMAP_1 <0)
Epi <- subset(x = Epi, subset = UMAP_2 <5)
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
cds <- order_cells(cds)

g <- plot_cells(cds, color_cells_by="cell.type",cell_size = 1,
                show_trajectory_graph = TRUE,trajectory_graph_segment_size = 0.75,
                label_cell_groups = F,
                label_groups_by_cluster = F, 
                label_branch_points = F)+scale_colour_manual(values = ExtractMetaColor(object))


jpeg(paste0(path,"pseudotime~.jpeg"), units="in", width=10, height=7,res=600)
print(g)+ ggtitle("Pseudotime trajectory from Seurat")+ 
        theme(plot.title = element_text(size=20, hjust = 0.5,face="plain"))
dev.off()
