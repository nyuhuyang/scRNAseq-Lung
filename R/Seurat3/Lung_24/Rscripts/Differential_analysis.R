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

(load(file = paste0("data/Lung_24_20191105.Rda")))
(load(file = paste0("data/Lung_24_20190918.Rda")))

# Differential analysis
Idents(object) = "integrated_snn_res.0.8"
Idents(object) = "RNA_snn_res.0.8"
DefaultAssay(object)  = "SCT"
object %<>% sortIdent(numeric = T)
Lung_markers <- FindAllMarkers.UMI(object, logfc.threshold = 0.5, only.pos = T,
                                   test.use = "MAST")
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
write.csv(Lung_markers,paste0(path,"Lung_24_Orignal-FC0.5_markers_res=0.8.csv"))
#Top_n = 5
#top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
#object %<>% ScaleData(features=unique(c(as.character(top$gene))))

#DoHeatmap.1(object, marker_df = Lung_markers, Top_n = 5, do.print=T,
#            angle = 0,group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
#            assay = "SCT",
#            label=T, cex.row=5, legend.size = 5,width=10, height=7,unique.name = T,
#            title = paste("Top 5 markers of each clusters in 24 sampels"))