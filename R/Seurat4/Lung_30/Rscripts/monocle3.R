# conda activate r4.1.1
#library(Seurat)
#library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(data.table)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# Need 16GB ?
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))
set.seed(args*1234)

step = c("monocle3","figures")[1]

if(step == "monocle3"){
    object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
    
    object@meta.data %<>% cbind(object[["umap"]]@cell.embeddings)
    Region = ifelse(args <= 10, "distal","COPD") 
    object <- subset(object, subset = Family %in% c("ASE","AT")
                     & Regions == Region
                     & UMAP_2 < -2)
    
    #============== re-run UMAP ==============================
    # Building trajectories with Monocle 3
    cds <- new_cell_data_set(object[["SCT"]]@data,
                             cell_metadata = object@meta.data)
    rm(object);GC()
    
    meta.data = pData(cds)
    meta.data = meta.data[!duplicated(meta.data$Cell_subtype),
                          c("Cell_subtype","Cell_subtype.colors")]
    meta.data = meta.data[order(meta.data$Cell_subtype),]
    # Run PCA
    cds <- preprocess_cds(cds, method = "PCA", num_dim = (100-args))
    cds <- reduce_dimension(cds, reduction_method = "UMAP",
                            preprocess_method = "PCA", init = "random",
                            umap.fast_sgd = FALSE)
    cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
    cds <- learn_graph(cds, use_partition = FALSE)
    
    get_earliest_principal_node <- function(cds, cell.type="BC"){
        cell_ids <- which(colData(cds)[, "Cell_subtype"] == cell.type)
        
        closest_vertex <-
            cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
        root_pr_nodes <-
            igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                      (which.max(table(closest_vertex[cell_ids,]))))]
        
        root_pr_nodes
    }
    cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
    
    #=============== prepare figures ===================
    # plot trajectories colored by pseudotime
    num = args
    if(num < 10) num = paste0("0",num)
    save.path = paste0(path,"rerun UMAP_",Region,"_",num,"/")
    if(!dir.exists(save.path))dir.create(save.path, recursive = T)
    
    for(label in c("subtype_raw","subtype","subtype_label","subtype_legend")){
        g <- plot_cells(
            cds = cds,
            label_cell_groups = switch(label,
                                       "subtype_raw" = TRUE,
                                       "subtype" = TRUE,
                                       "subtype_label" = TRUE,
                                       "subtype_legend" = FALSE),
            color_cells_by = "Cell_subtype",
            show_trajectory_graph = switch(label,
                                           "subtype_raw" = FALSE,
                                           TRUE),
            group_label_size = switch(label,
                                      "subtype_raw" = 5,
                                      "subtype_label" = 5,
                                      0)) +
                scale_color_manual(values=meta.data[,"Cell_subtype.colors"])
        jpeg(paste0(save.path,Region,"_Cell_",label,".jpeg"), units="in", width=7, height=7,res=600)
        print(g)
        dev.off()
    }
    
    for(label in c("pseudotime","pseudotime_label")){
        g1 <-     plot_cells(
                    cds = cds,
                    color_cells_by = "pseudotime",
                    show_trajectory_graph = switch(label,
                                                   "pseudotime" = FALSE,
                                                   "pseudotime_label" = TRUE
                                                   )
                    )
        jpeg(paste0(save.path,Region,"_",label,".jpeg"), units="in", width=7, height=7,res=600)
        print(g1)
        dev.off()
    }
        

    
    # save coordinates
    pData(cds)$pseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
    
    fwrite(as.data.frame(pData(cds)[c("barcode","pseudotime","UMAP_1","UMAP_2")]),
           file = paste0(save.path,Region,"_pseudotime.csv"))
    fwrite(as.data.frame(cds@preprocess_aux@listData$gene_loadings),
           file = paste0(save.path,Region,"_PCA.csv"),row.names = TRUE)
    saveRDS(cds,file = paste0(save.path,"cds_",Region,"_",num,".rds"))
}

