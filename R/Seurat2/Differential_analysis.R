########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kableExtra)
library(SingleR)
library(reshape2)
library(MAST)
source("../R/Seurat_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 4.1 load data & Rename ident, Compare DE across all major cell types=========================
# 4.1.1 load data ==============
(load(file = "./data/Lung_MNN_9_20181101.Rda"))
table(object@meta.data$orig.ident)
table(object@ident)

TSNEPlot.1(object = object,do.label = T, group.by = "ident",
            do.return = TRUE, no.legend = T,colors.use = ExtractMetaColor(object),
            pt.size = 1,label.size = 5,label.repel = T)+
            ggtitle("Manually label cell types")+
            theme(text = element_text(size=20),
            plot.title = element_text(hjust = 0.5, face = "bold"))

# 4.1.2 Compare DE across all cell types ==============
object %<>% SetAllIdent(id = "res.0.8")
markers.res.0.8 <- FindAllMarkers.UMI(object,test.use = "MAST",
                                  logfc.threshold = 0.1,return.thresh = 0.05)
write.csv(markers.res.0.8,paste0(path,"markers.res.0.8.csv"))
write.csv(as.matrix(object@data),paste0(path,"All_expression.csv"))


Epi %<>% SetAllIdent(id = "res.0.8")
markers.Epi <- FindAllMarkers.UMI(Epi,test.use = "MAST",
                                  logfc.threshold = 0.1,return.thresh = 0.05)
write.csv(markers.Epi,paste0(path,"markers.Epi.csv"))
write.csv(as.matrix(Epi@data),paste0(path,"Epi_expression.csv"))

# 4.1.3 Compare DE across all samples ==============

markers_list <- lapply(Split.object, function(sub_object) {
        FindAllMarkers.UMI(sub_object,test.use = "MAST",
                           logfc.threshold = 0.25,return.thresh = 0.05)
        })
names(markers_list) = names(Split.object)
lapply(names(markers_list), function(s) {
        write.csv(markers_list[s],paste0(path,"markers_",s,".csv"))
})

for(i in 1:length(Split.object)){
        write.csv(as.matrix(Split.object[[i]]@data),
                  paste0(path,unique(Split.object[[i]]@meta.data$orig.ident),
                         "_expression.csv"))  
}

 for(i in 1:3){
        DoHeatmap.1(Split.object[[i]], markers_list[[i]], Top_n = 10,
                    ident.use = unique(Split.object[[i]]@meta.data$conditions),group.label.rot = F,cex.row = 4,remove.key =F,
                    title.size = 18,do.print = T)
}
object@scale.data = readRDS("data/Lung.scale.data_Harmony_3_20190109.rds")
DoHeatmap.1(object, markers.DPT, Top_n = 10,
            ident.use = paste(unique(object@meta.data$conditions),collapse = "_"),
            group.label.rot = F,cex.row = 3,remove.key =F,
            title.size = 18,do.print = T)

