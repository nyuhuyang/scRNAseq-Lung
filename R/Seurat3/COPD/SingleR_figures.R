library(Seurat)
library(magrittr)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggthemes)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file = "data/Lung_16_distal_20191022.Rda"))
(load(file="output/singlerT_Lung_16_distal_20191023.Rda"))
singler@nrows == ncol(object)
table(rownames(singler) == colnames(object))
##############################
# add singleR label to Seurat
###############################

singlerDF = data.frame("labels" = as.character(singler$pruned.labels),
                       row.names = singler@rownames, stringsAsFactors =F)
sapply(singlerDF,class)
singlerDF$labels[is.na(singlerDF$labels)] = "Unknown"
head(singlerDF)
table(object$labels) %>% kable() %>% kable_styling()
object <- AddMetaData(object = object,metadata = singlerDF)
DefaultAssay(object) <- 'RNA'
object %<>% FindClusters(resolution = 1)
object %<>% FindClusters(resolution = 0.7)
object@meta.data[object$RNA_snn_res.1 == 4,"labels"] = "Fibroblasts:TBX4"
object@meta.data[object$RNA_snn_res.1 == 17,"labels"] = "Fibroblasts:IGF1"
object@meta.data[object$RNA_snn_res.1 == 20,"labels"] = "Fibroblasts:CCL19"
object@meta.data[object$RNA_snn_res.1 == 2,"labels"] = "T cells:TRM"
object@meta.data[object$RNA_snn_res.1 == 7,"labels"] = "T cells:CD4+"
object@meta.data[object$RNA_snn_res.1 == 8,"labels"] = "Basal cells"
object@meta.data[object$RNA_snn_res.1 == 27,"labels"] = "Intermediate cells"
object@meta.data[object$RNA_snn_res.1 == 21,"labels"] = "Neutrophils"
object@meta.data[object$RNA_snn_res.1 %in% c(10,18),"labels"] = "Alveolar type 2"

object@meta.data[object$RNA_snn_res.1 %in% c(5,31) & 
                         object$labels != "Secretory cells:Distal" ,"labels"] = "Secretory cells"
object@meta.data[object$RNA_snn_res.1 == 14 & object$labels !="Dendritic cells" ,
                 "labels"] = "Macrophages"

object@meta.data[object$RNA_snn_res.0.7 == 0,"labels"] = "Ciliated cells:1"
object@meta.data[object$RNA_snn_res.0.7 == 9,"labels"] = "Ciliated cells:2"
object@meta.data[object$RNA_snn_res.0.7 == 20,"labels"] = "Ciliated cells:3"
object@meta.data[object$RNA_snn_res.0.7 == 26,"labels"] = "Pre-ciliated cells"

object@meta.data = cbind(object@meta.data, object@reductions$umap@cell.embeddings)
object@meta.data[(object$UMAP_1 < 0 & object$UMAP_2 >10 & object$RNA_snn_res.1 == 35),
                 "labels"] = "Alveolar type 2:Proliferating"

object@meta.data[(object$UMAP_2 > -2 &  grepl(".*gland cells",object$labels)),
                 "labels"] = "Secretory cells"
object@meta.data[(object$UMAP_2 < -2 &  grepl(".*gland cells",object$labels)),
                 "labels"] = "Unknown"

Idents(object) = "labels"
SK <- subset(object, idents= "T cells:7SK.2+")

UMAPPlot.1(SK, group.by = "RNA_snn_res.0.7",label = T)
Idents(SK) = "RNA_snn_res.0.7"
SK %<>%  RenameIdents("1"="Monocytes",
                      "3"="NK cells",
                      "4"="Endothelial cells:HEVs",
                      "5"="Endothelial cells:Capillary",
                      "6"="T cells:7SK.2+",
                      "8"="Smooth muscle:Vascular",
                      "13"="B cells",
                      "14"="Alveolar macrophages",
                      "16"="Mast cells",
                      "17"="Neutrophils",
                      "18"="Smooth muscle:Airway",
                      "22"="Endothelial cells:Arterial",
                      "24"="Endothelial cells:Lymphatic",
                      "27"="B cells:Plasma",
                      "28"="Endothelial cells:Proliferating",
                      "29"="Neurons",
                      "30"="Alveolar type 1",
                      "31"="Ionocytes/NEC",
                      "32"="T cells:7SK.2+")
SK@meta.data$labels = as.character(Idents(SK))
object@meta.data[colnames(SK),"labels"] = SK@meta.data$labels
object@meta.data = object@meta.data[,-grep("UMAP",colnames(object@meta.data))]
Idents(object) <- "labels"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "labels", colors = Singler.colors)
UMAPPlot.1(object, cols = ExtractMetaColor(object),label = T, label.repel = T,pt.size = 0.3,
           label.size = 3, repel = T,no.legend = T,do.print = T,
           do.return = F,alpha = 0.9,
           title = "Minor cell types",unique.name = "conditions")

save(object,file="data/Lung_16_distal_20191019.Rda")

(load(file = "data/Lung_16_distal_20191019.Rda"))

meta.data = object@meta.data[,c("labels","labels.colors")]
meta.data$Barcode = rownames(meta.data)
