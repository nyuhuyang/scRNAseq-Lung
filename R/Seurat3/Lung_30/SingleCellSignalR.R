#conda activate r4.0
invisible(lapply(c("dplyr","SingleCellSignalR","cowplot","openxlsx",
                   "magrittr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data
object = readRDS(file = "data/Lung_30_20200710.rds")
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
object %<>% AddMetaColor(label= "annotations3", colors = Singler.colors)
Idents(object) = "conditions"

conditions <- c("proximal","distal","terminal","COPD")[-1]
for (con in  conditions){
    save.path = paste0(path,"SingleCellSignalR/",con,"/")
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
    setwd(save.path)
    print(pracma::pwd())
    sub_object <- subset(object, idents = con)
    #Idents(sub_object) = "orig.ident"
    #sub_object %<>% subset(idents = "UNC-48-P")
    
    df<- data.frame("cluster" = as.numeric(as.factor(sub_object$annotations3)),
                    "annotations3" = sub_object$annotations3,
                    "cluster_num" = paste("cluster",as.numeric(as.factor(sub_object$annotations3))))
    data = data.frame(sub_object[["SCT"]]@data)
    all.genes <- rownames(data)
    
    # Ligand/Receptor analysis using SingleCellSignalR
    signal = cell_signaling(data=data,genes=all.genes,cluster=df$cluster)
    
    # re-name cluster name to cell types
    df_short = df[!duplicated(df$cluster),]
    
    signal_annotations3 = signal
    signal_annotations3 %<>% lapply(function(x){
        names(x) %<>% plyr::mapvalues(from = df_short$cluster_num, to = df_short$annotations3,
                                      warn_missing = FALSE)
        x
    })
    name_pairs =names(signal_annotations3)
    df_name_pairs = stringr::str_split(name_pairs, pattern = "-")
    
    df_name_pairs %<>% sapply(function(x){
        x %<>% plyr::mapvalues(from = df_short$cluster_num, to = df_short$annotations3,
                                      warn_missing = FALSE)
        paste0(x[1],"-",x[2])
    })
    names(signal_annotations3) = df_name_pairs
    
    signal_annotations3_short <- lapply(signal_annotations3, function(x){
        x[x$LRscore > 0.968,]
    })
    signal_annotations3_short = signal_annotations3_short[lapply(signal_annotations3_short,nrow)>0]
    length(signal_annotations3_short)
    # Visualization
    jpeg(paste0("visualize_interactions_",con,".jpeg"), 
         units="in", width=5, height=5,res=600)
    visualize_interactions(signal_annotations3_short)
    dev.off()
    jpeg(paste0("visualize_interactions_show.in.jpeg"), 
         units="in", width=5, height=5,res=600)
    visualize_interactions(signal_annotations3, show.in = c(1,4))
    dev.off()
    
    # save txt
    if(!dir.exists("cell-signaling_annotation/")) dir.create("cell-signaling_annotation/", recursive = T)
    for(i in 1:length(signal_annotations3)){
        write.table(signal_annotations3[[i]],file = paste0("cell-signaling_annotation/LR_interactions_",
                                                           names(signal_annotations3)[i],
                                                           "-paracrine.txt"))
        Progress(i,length(signal_annotations3))
    }
    setwd("../../../..")
}


intra = intra_network("S1PR1",data,all.genes,cluster,"cluster 3",signal = signal)
