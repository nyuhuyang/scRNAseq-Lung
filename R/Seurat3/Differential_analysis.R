########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
library(MAST)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file = "data/Lung_harmony_12_20190614.Rda"))
Idents(object) <- "RNA_snn_res.1.2"
Lung_markers <- FindAllMarkers.UMI(object, logfc.threshold = 0.1,
                                   only.pos = T)

write.csv(Lung_markers,paste0(path,"Lung_12_markers.csv"))
#Lung.markers.csv =read.csv(file = paste0(path,"Lung.markers.csv"),
#                                row.names = 1, stringsAsFactors=F)
TSNEPlot.1(object,label = F, repel = F, no.legend = F,pt.size = 1,
           cols = ExtractMetaColor(object),do.return = T,do.print = F,
           title = "All clusters in in CU12-D, CU12-D-repeat, and CU12-T")

object %<>% ScaleData(features=unique(Lung_markers$gene))
DoHeatmap.1(object, marker_df = Lung_markers, Top_n = 5, do.print=T, angle = 0,
            group.bar = T, title.size = 13, no.legend = F,size=0,hjust = 0.5,
            label=F, cex.row=5, legend.size = 5,width=10, height=7,
            title = "Top 5 markers in each clusters")

#split by samples================
Idents(object) <- 'orig.ident'
table(Idents(object))

#CU12-D
CU12_D <- subset(object, ident = "CU12-D")
Idents(CU12_D) <-"singler1sub"

TSNEPlot.1(CU12_D,label = F, repel = F, no.legend = F,pt.size = 1,
           cols = ExtractMetaColor(CU12_D),do.return = F,do.print = T,
           title = "All clusters in CU12-D")

CU12_D.markers <- FindAllMarkers.UMI(object = CU12_D, logfc.threshold = 1,
                                     only.pos = T,test.use = "MAST")
write.csv(CU12_D.markers,paste0(path,"CU12_D.markers.csv"))
table(CU12_D.markers$cluster)
DoHeatmap.1(CU12_D, CU12_D.markers, Top_n = 10, do.print=T, angle = 0,
            group.bar = T, title.size = 13, no.legend = F,size=0,hjust = 0.5,
            label=T, cex.row=4, legend.size = 5,width=10, height=7,
            title = "Top 10 markers in all clusters in CU12-D")

#CU12-D-repeat
CU12_D_r <- subset(object, ident = "CU12-D-Repeat")
Idents(CU12_D_r) <-"singler1sub"

TSNEPlot.1(CU12_D_r,label = F, repel = F, no.legend = F,pt.size = 1,
           cols = ExtractMetaColor(CU12_D_r),do.return = F,do.print = T,
           title = "All clusters in CU12-D-Repeat")
CU12_D_r.markers <- FindAllMarkers.UMI(object = CU12_D_r, logfc.threshold = 1,
                                     only.pos = T,test.use = "MAST")
write.csv(CU12_D_r.markers,paste0(path,"CU12_D_r.markers.csv"))

DoHeatmap.1(CU12_D_r, CU12_D_r.markers, Top_n = 10, do.print=T, angle = 0,
            group.bar = T, title.size = 13, no.legend = F,size=0,hjust = 0.5,
            label=T, cex.row=4, legend.size = 5,width=10, height=7,
            title = "Top 3 markers in all clusters in CU12-D-Repeat")

#CU12-T
CU12_T <- subset(object, ident = "CU12-T")
Idents(CU12_T) <-"singler1sub"

TSNEPlot.1(CU12_T,label = F, repel = F, no.legend = F,pt.size = 1,
           cols = ExtractMetaColor(CU12_T),do.return = F,do.print = T,
           title = "All clusters in CU12-T")
CU12_T.markers <- FindAllMarkers.UMI(object = CU12_T, logfc.threshold = 1,
                                       only.pos = T,test.use = "MAST")
write.csv(CU12_T.markers,paste0(path,"CU12_T.markers.csv"))

table(CU12_T.markers$cluster)
DoHeatmap.1(CU12_T, CU12_T.markers, Top_n = 10, do.print=T, angle = 0,
            group.bar = T, title.size = 13, no.legend = F,size=0,hjust = 0.5,
            label=T, cex.row=4, legend.size = 5,width=10, height=7,
            title = "All clusters in CU12-T")
