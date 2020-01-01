####################################
invisible(lapply(c("Seurat","dplyr","monocle3","cowplot","magrittr",
                   "tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# read file
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
# Load Seurat

(load(file="data/Epi_24_All_20191223.Rda"))
table(object$orig.ident)
object$cell_types <- plyr::mapvalues(object$cell.types,
                                       from = df_cell_types$`Cell types`,
                                       to = df_cell_types$Abbreviation)
object$cell_types.colors = object$cell.types.colors

npcs =100
object <- RunPCA(object, features = VariableFeatures(object),verbose =F,npcs = npcs)

object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
#object %<>% RunTSNE(reduction = "pca", dims = 1:npcs, check_duplicates = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

object@assays$RNA@scale.data = matrix(0,0,0)

Idents(object) = "cell_types"
object %<>% sortIdent()
lapply(c(T,F),function(x) {
        UMAPPlot.1(object, group.by="cell_types",pt.size = 1,label = x,
                     label.repel = T,alpha = 0.9,cols = ExtractMetaColor(object),
                     no.legend = x,label.size = 4, repel = T, title = paste("All Epithelial cell types"),
                     do.return = F,do.print = T,unique.name = "conditions")}
        )

#Construct monocle cds
cds <- new_cell_data_set(object@assays$SCT@counts,
                         cell_metadata = object@meta.data,
                         gene_metadata = data.frame("gene_short_name" = rownames(object),
                                                      row.names = rownames(object)))
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = npcs)
## Step 2: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)
reducedDims(cds)$UMAP <- object@reductions$umap@cell.embeddings

## Step 3: Cluster the cells
cds <- cluster_cells(cds)

## Step 4: Learn a graph
cds <- learn_graph(cds)

cell_meata <- colData(cds)
if(table(rownames(cell_meata) == colnames(object))) {
        cell_meata$cell_types = object$cell_types
        cell_meata$cell_types.colors = object$cell_types.colors
}
colData(cds) = cell_meata

#Idents(object) = "cell_types"
#object %<>% sortIdent
#g <- plot_cells(cds, color_cells_by="cell_types",cell_size = 1,
#                label_cell_groups = F,
#                label_groups_by_cluster = F, 
#                label_branch_points = F)+scale_colour_manual(values = ExtractMetaColor(object))
#jpeg(paste0(path,"Plot_Epi_cell.type.jpeg"), units="in", width=10, height=7,res=600)
#print(g)+ ggtitle(paste("Trajectory of",con,"epithelial cell types"))+ 
#        theme(plot.title = element_text(size=20, hjust = 0.5,face="plain"))
#dev.off()

## Step 5: Order cells
#roots <- rownames(cell_meata)[cell_meata$cell_types %in% c("BC","MEC","Sec-d2")]
#cds <- order_cells(cds, root_cells=roots)
#g <- plot_cells(cds, color_cells_by="pseudotime",cell_size = 1,
#                label_cell_groups = F,
#                label_groups_by_cluster = F, 
#                label_branch_points = F,
#                label_roots = F, label_leaves = T)
#jpeg(paste0(path,"Plot_Epi_",con,"_pseudotime~.jpeg"), units="in", width=10, height=7,res=600)
#print(g)+ ggtitle(paste("Pseudotime trajectory of",con,"epithelial cells"))+ 
#        theme(plot.title = element_text(size=20, hjust = 0.5,face="plain"))
#dev.off()
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin=c("BC","MEC","Sec-d2")){
        root_pr_nodes <- NULL
        for(k in seq_along(time_bin)){
                cell_ids <- which(colData(cds)[, "cell_types"] %in% time_bin[k])
                closest_vertex <-
                        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
                closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
                top <- which.max(table(closest_vertex[cell_ids,]))
                
                root_pr_nodes <- c(root_pr_nodes,
                                  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(top))])
        }

        root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds,time_bin=c("BC")))
save(object,cds, file = paste0("data/Epi_24_All_20191223.Rda"))

# subset by region

conditions = c("All","proximal","distal","terminal")
Idents(object) = "cell.types"
for (i in seq_along(conditions) ){
        con = conditions[i]
        if (con == "proximal") {
                sub_object <- subset(object, idents = c("Alveolar type 1 cells","Alveolar type 2 cells:A",
                                              "Alveolar type 2 cells:B","Alveolar type 2 cells:C",
                                              "Secretory cells:Distal","Secretory cells:Distal:2"), 
                                   invert = TRUE)
        }
        if (con == "distal") {
                sub_object <- subset(object, idents = c("Myoepithelial cells","Squamous",
                                              "Submucosal gland:Mucous cells",
                                              "Submucosal gland:Serous cells"), 
                                   invert = TRUE)
        }
        if (con == "terminal") {
                sub_object <- subset(object, idents = c("Mucus-producing cells","Myoepithelial cells",
                                              "Squamous","Submucosal gland:Mucous cells",
                                              "Submucosal gland:Serous cells"),
                                   invert = TRUE)
        }
        if (con == "All") sub_object <- object
        Idents(sub_object) = "conditions"
        if(con != "All") sub_object %<>% subset(idents = con)
        sub_cds <- cds[,colnames(sub_object)]
        g <- plot_cells(sub_cds, color_cells_by="pseudotime",cell_size = 1,
                        label_cell_groups = F,
                        label_groups_by_cluster = F, 
                        label_branch_points = F)+ ggtitle(paste("Pseudotime trajectory of",con,"epithelial cells"))+ 
                theme(plot.title = element_text(size=20, hjust = 0.5,face="plain"))
        jpeg(paste0(path,"Plot_Epi_",con,"_pseudotime.jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        
        Idents(sub_object) = "cell_types"
        sub_object %<>% sortIdent()
        lapply(c(T,F),function(x) {
                UMAPPlot.1(sub_object, group.by="cell_types",pt.size = 1,label = x,
                           label.repel = T,alpha = 0.9,cols = ExtractMetaColor(object),
                           no.legend = x,label.size = 4, repel = T, title = paste("All Epithelial cell types"),
                           do.return = F,do.print = T,unique.name = "conditions")}
        )
        
        Progress(i, length(conditions))
}



