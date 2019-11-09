.p .,                                                             ########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(data.table)
library(kableExtra)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
(load(file="data/Lung_16_distal_20191022.Rda"))
distal_label <- object@meta.data[,c("labels","Doublets")]
colnames(distal_label)[1] = "cell.types"
distal_label$barcode = rownames(distal_label)
(load(file="data/Lung_24proximal_20191004.Rda"))
proximal_label <- object@meta.data[,c("cell.types","Doublets")]
proximal_label$barcode = rownames(proximal_label)
(load(file="data/Lung_24terminal_20191004.Rda"))
terminal_label <- object@meta.data[,c("cell.types","Doublets")]
terminal_label$barcode = rownames(terminal_label)
df_Labels <- do.call(rbind,list(distal_label,proximal_label,terminal_label))

(load(file = "data/Lung_24_20191105.Rda"))
meta.data = object@meta.data %>% mutate(barcode = rownames(.))
meta.data = left_join(meta.data,df_Labels,by = "barcode")
rownames(meta.data) = meta.data$barcode
meta.data = meta.data[colnames(object),]
table(meta.data$cell.types)
meta.data$cell.types %<>% gsub("arterial","Arterial",.)
meta.data$cell.types %<>% gsub("^Endothelial cells:HEV$","Endothelial cells:HEVs",.)
meta.data$cell.types %<>% gsub("^Endothelial cells:smooth muscel$","Endothelial cells:Smooth muscle",.)
meta.data$cell.types %<>% gsub("^Endothelial:Smooth muscle$","Endothelial cells:Smooth muscle",.)
meta.data$cell.types %<>% gsub("^Hybrid$","Hybrid cells",.)
meta.data$cell.types %<>% gsub("^Macrophage:Prolifereating$","Macrophages:Proliferating",.)
meta.data$cell.types %<>% gsub("^Neuro endocrine$","Neuroendocrine",.)
meta.data$cell.types %<>% gsub("^Neutrophil$","Neutrophils",.)
meta.data$cell.types %<>% gsub("^Proliferating basal cells$","Basal cells:Proliferating",.)
meta.data$cell.types %<>% gsub("^T cells:7SK.2$","T cells:7SK.2+",.)

object@meta.data = meta.data
Idents(object) = "cell.types"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "cell.types", colors = Singler.colors)
UMAPPlot.1(object, cols = ExtractMetaColor(object),label = T, label.repel = T,pt.size = 1,
           label.size = 3, repel = T,no.legend = T,do.print = T,
           title = "Cell types")
TSNEPlot.1(object, cols = ExtractMetaColor(object),label = T, label.repel = T,pt.size = 1,
           label.size = 3, repel = T,no.legend = T,do.print = T,
           title = "Cell types")
