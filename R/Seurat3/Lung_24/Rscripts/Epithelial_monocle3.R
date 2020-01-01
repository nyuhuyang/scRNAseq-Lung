####################################
invisible(lapply(c("Seurat","dplyr","monocle3","cowplot","magrittr",
                   "tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
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

# conditions
conditions = c("proximal","distal","terminal","All")
(con <- conditions[args])

# read file
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
# Load Seurat

(load(file = paste0("data/Epi_24_",con,"_20191223.Rda")))
object$cell_types <- plyr::mapvalues(object$cell.types,
                                     from = df_cell_types$`Cell types`,
                                     to = df_cell_types$Abbreviation)
object$cell_types.colors = object$cell.types.colors
Idents(object) = "cell_types"

UMAPPlot.1(object, group.by="cell_types",pt.size = 1,label = F,
           label.repel = T,alpha = 0.9,cols = ExtractMetaColor(object),
           no.legend = F,label.size = 4, repel = T, title = paste(con,"cell types"),
           do.return = F,do.print = F,unique.name = "conditions")

#Construct monocle cds
#cds <- new_cell_data_set(object@assays$SCT@counts,
#                         cell_metadata = object@meta.data,
#                         gene_metadata = data.frame("gene_short_name" = rownames(object),
row.names = rownames(object)))
## Step 1: Normalize and pre-process the data
#cds <- preprocess_cds(cds, num_dim = 50)
## Step 2: Reduce the dimensions using UMAP
#cds <- reduce_dimension(cds)
## Step 3: Cluster the cells
#cds <- cluster_cells(cds)

## Step 4: Learn a graph
#cds <- learn_graph(cds)

df <- colData(cds)
if(table(rownames(df) == colnames(object))) {
        df$cell_types = object$cell_types
        df$cell_types.colors = object$cell_types.colors
}
colData(cds) = df


Idents(object) = "cell_types"
object %<>% sortIdent
g <- plot_cells(cds, color_cells_by="cell_types",cell_size = 1,
                label_cell_groups = F,
                label_groups_by_cluster = F, 
                label_branch_points = F)+scale_colour_manual(values = ExtractMetaColor(object))
jpeg(paste0(path,"Plot_Epi_",con,"_cell.type.jpeg"), units="in", width=10, height=7,res=600)
print(g)+ ggtitle(paste("Trajectory of",con,"epithelial cell types"))+ 
        theme(plot.title = element_text(size=20, hjust = 0.5,face="plain"))
dev.off()

## Step 5: Order cells
cds <- order_cells(cds)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="BC"){
        cell_ids <- which(colData(cds)[, "cell_types"] == time_bin)
        
        closest_vertex <-
                cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
        root_pr_nodes <-
                igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                          (which.max(table(closest_vertex[cell_ids,]))))]
        
        root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
g <- plot_cells(cds, color_cells_by="pseudotime",cell_size = 1,
                label_cell_groups = T,
                label_groups_by_cluster = F, 
                label_branch_points = F)
jpeg(paste0(path,"Plot_Epi_",con,"_pseudotime.jpeg"), units="in", width=10, height=7,res=600)
print(g)+ ggtitle(paste("Pseudotime trajectory of",con,"epithelial cells"))+ 
        theme(plot.title = element_text(size=20, hjust = 0.5,face="plain"))
dev.off()

save(object,cds, file = paste0("data/Epi_24_",con,"_20191223.Rda"))
