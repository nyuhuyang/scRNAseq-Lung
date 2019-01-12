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
source("./R/Seurat_functions.R")

path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
dir.create(path, recursive = T)
# 4.1 load data & Rename ident, Compare DE across all major cell types=========================
# 4.1.1 load data ==============
lname1 = load(file = "./data/object_20180825.Rda");lname1object@meta.data$orig.ident = gsub("PND18pre","PND18",object@meta.data$orig.ident)
table(object@meta.data$orig.ident)
table(object@ident)

TSNEPlot.1(object = object,do.label = T, group.by = "ident",
            do.return = TRUE, no.legend = T,colors.use = singler.colors,
            pt.size = 1,label.size = 5,label.repel = T)+
            ggtitle("Manually label cell types")+
            theme(text = element_text(size=20),
            plot.title = element_text(hjust = 0.5, face = "bold"))
Split.object[[1]] <- SetAllIdent(Split.object[[1]], id = "res.0.6")
Split.object[[2]] <- SetAllIdent(Split.object[[2]], id = "res.0.6")
Split.object[[3]] <- SetAllIdent(Split.object[[3]], id = "res.0.6")

# 4.1.2 Compare DE across all cell types ==============
markers_list <- lapply(Split.object, function(sub_object) {
        FindAllMarkers.UMI(sub_object,test.use = "MAST",
                           logfc.threshold = 2,return.thresh = 0.05)
        })

markers.D <- FindAllMarkers.UMI(Split.object[[1]],test.use = "MAST",
                                    logfc.threshold = 0.25,return.thresh = 0.05)
markers.P <- FindAllMarkers.UMI(Split.object[[2]],test.use = "MAST",
                                    logfc.threshold = 0.25,return.thresh = 0.05)
markers.T <- FindAllMarkers.UMI(Split.object[[3]],test.use = "MAST",
                                    logfc.threshold = 0.25,return.thresh = 0.05)
markers.DPT <- FindAllMarkers.UMI(object,test.use = "MAST",
                                      logfc.threshold = 0.25,return.thresh = 0.05)

write.csv(markers.D,paste0(path,"markers.D.csv"))
write.csv(markers.P,paste0(path,"markers.P.csv"))
write.csv(markers.T,paste0(path,"markers.T.csv"))
write.csv(markers.DPT,paste0(path,"markers.DPT.csv"))

write.csv(as.matrix(Split.object[[1]]@data),paste0(path,"D_expression.csv"))
write.csv(as.matrix(Split.object[[2]]@data),paste0(path,"P_expression.csv"))
write.csv(as.matrix(Split.object[[3]]@data),paste0(path,"T_expression.csv"))
write.csv(as.matrix(object@data),paste0(path,"DPT_expression.csv"))

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

