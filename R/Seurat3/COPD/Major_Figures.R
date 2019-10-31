########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(magrittr)
library(DoubletFinder)
library(kableExtra)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# clustering dendrogram (Spearman or Pearson) for all cell types/clusters (16 D+COPD) 
# based on top (50-100?) genes in each clust
(load(file = "data/Lung_16_distal_20191022.Rda"))
Lung_markers = read.csv(file = "Yang/distal_COPD/DE analysis/Lung_16_RNA_snn_res.0.8_markers.csv",
                        row.names = 1,stringsAsFactors = F)
Idents(object) = "RNA_snn_res.0.8"
table(Idents(object))
Top_n = 50
top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
object %<>% ScaleData(features=unique(c(as.character(top$gene))))
features = c(as.character(top$gene))
featuresNum <- make.unique(features, sep = ".")
object %<>% MakeUniqueGenes(features = features)

DoHeatmap.1(object, features = featuresNum, Top_n = Top_n, do.print=T,
            angle = 0,group.bar = T, title.size = 13, no.legend = T,size=4,hjust = 0.5,
            assay = "RNA",
            label=T, cex.row= 0.1, legend.size = 5,width=10, height=7,unique.name = T,
            title = paste("Top",Top_n,"markers of each cluster in 16 D+COPD sampels"))
