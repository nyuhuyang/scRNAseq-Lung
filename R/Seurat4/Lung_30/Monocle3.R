# conda activate r4.1.1
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(data.table)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#=========== sample 30 ==================
object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
object@meta.data %<>% cbind(object[["umap"]]@cell.embeddings)
D <- subset(object, subset = Family %in% c("ASE","AT")
            &  Regions == "distal" 
            & UMAP_2 < -2)
COPD <- subset(object, subset = Family %in% c("ASE","AT")
               &  Regions == "COPD"
               & UMAP_2 < -2)
#============== re-run UMAP ==============================
# Building trajectories with Monocle 3
for(i in 1:10){
    print(paste0("---------",i,"--------------"))
    set.seed(123*i)
    D.cds <- as.cell_data_set(D)
    meta.data1 = pData(D.cds)
    meta.data1 = meta.data1[!duplicated(meta.data1$Cell_subtype),
                            c("Cell_subtype","Cell_subtype.colors")]
    meta.data1 = meta.data1[order(meta.data1$Cell_subtype),]
    D.cds %<>% reduce_dimension(umap.fast_sgd = TRUE)
    D.cds <- cluster_cells(cds = D.cds, reduction_method = "UMAP")
    D.cds <- learn_graph(D.cds, use_partition = FALSE)
    D.cds <- order_cells(D.cds)
    
    COPD.cds <- as.cell_data_set(COPD)
    meta.data2 = pData(COPD.cds)
    meta.data2 = meta.data2[!duplicated(meta.data2$Cell_subtype),
                            c("Cell_subtype","Cell_subtype.colors")]
    meta.data2 = meta.data2[order(meta.data1$Cell_subtype),]
    COPD.cds %<>% reduce_dimension(umap.fast_sgd = TRUE)
    COPD.cds <- cluster_cells(cds = COPD.cds, reduction_method = "UMAP")
    COPD.cds <- learn_graph(COPD.cds, use_partition = FALSE)
    COPD.cds <- order_cells(COPD.cds)
    
    #=============== prepare figures ===================
    # plot trajectories colored by pseudotime
    num = i
    if(num < 10) num = paste0("0",num)
    save.path = paste0(path,"rerun UMAP_",num,"/")
    if(!dir.exists(save.path))dir.create(save.path, recursive = T)
    
    for(label in c("subtype","subtype_label","subtype_legend")){
        jpeg(paste0(save.path,"distal_Cell_",label,".jpeg"), units="in", width=7, height=7,res=600)
        print(plot_cells(
            cds = D.cds,
            label_cell_groups = switch(label,
                                       "subtype" = TRUE,
                                       "subtype_label" = TRUE,
                                       "subtype_legend" = FALSE),
            color_cells_by = "Cell_subtype",
            show_trajectory_graph = TRUE,group_label_size = switch(label,
                                                                   "subtype_label" = 5,
                                                                   0)) +
                scale_color_manual(values=meta.data1[,"Cell_subtype.colors"]))
        dev.off()
    }
    
    
    jpeg(paste0(save.path,"distal_pseudotime.jpeg"), units="in", width=7, height=7,res=600)
    plot_cells(
        cds = D.cds,
        color_cells_by = "pseudotime",
        show_trajectory_graph = TRUE
    )
    dev.off()
    
    for(label in c("subtype","subtype_label","subtype_legend")){
        jpeg(paste0(save.path,"COPD_Cell_",label,".jpeg"), units="in", width=7, height=7,res=600)
        print(plot_cells(
            cds = COPD.cds,
            label_cell_groups = switch(label,
                                       "subtype" = TRUE,
                                       "subtype_label" = TRUE,
                                       "subtype_legend" = FALSE),
            color_cells_by = "Cell_subtype",
            show_trajectory_graph = TRUE,group_label_size = switch(label,
                                                                   "subtype_label" = 5,
                                                                   0)) +
                scale_color_manual(values=meta.data2[,"Cell_subtype.colors"]))
        dev.off()
    }
    
    
    jpeg(paste0(save.path,"COPD_pseudotime.jpeg"), units="in", width=7, height=7,res=600)
    plot_cells(
        cds = COPD.cds,
        color_cells_by = "pseudotime",
        show_trajectory_graph = TRUE
    )
    dev.off()
    
    # save coordinates
    pData(D.cds)$pseudotime <- D.cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
    pData(COPD.cds)$pseudotime <- COPD.cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
    
    fwrite(as.data.frame(pData(D.cds)[c("barcode","pseudotime","UMAP_1","UMAP_2")]),
           file = paste0(save.path,"distal_pseudotime.csv"))
    
    fwrite(as.data.frame(pData(COPD.cds)[c("barcode","pseudotime","UMAP_1","UMAP_2")]),
           file = paste0(save.path,"COPD_pseudotime.csv"))
}



