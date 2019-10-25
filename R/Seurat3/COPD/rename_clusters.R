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
meta.data1 = object@meta.data
meta.data1$Barcode = rownames(meta.data1)
meta.data2 = full_join(meta.data1,meta.data, by = "Barcode")
rownames(meta.data2) = meta.data2$Barcode
meta.data2[is.na(meta.data2$labels),"labels"] = "Unknown"

Idents(object) = "orig.ident"
VU19 <- subset(object, idents = "VU19-D")
UMAPPlot.1(VU19, group.by = "RNA_snn_res.0.8",label = T, do.print = T)

rename <- object@meta.data[!duplicated(object$RNA_snn_res.0.8),c("labels","RNA_snn_res.0.8")]
rownames(rename) = NULL
rename = rename[order(rename$RNA_snn_res.0.8),]
Idents(VU19) = "RNA_snn_res.0.8"
VU19 %<>% RenameIdents("0"="Fibroblasts:IGF1",
                       "1"="Ciliated cells:1",
                       "2"="Endothelial cells:HEVs",
                       "3"="NK cells",
                       "4"="T cells:TRM",
                       "5"="Monocytes",
                       "6"="Endothelial cells:Capillary",
                       "7"="Pericytes",
                       "8"="Secretory cells",
                       "9"="Ciliated cells:2",
                       "10"="Basal cells",
                       "11"="T cells:CD4+",
                       "12"="B cells",
                       "13"="Alveolar macrophages",
                       "14"="Ciliated cells:3",
                       "15"="Macrophages",
                       "16"="Alveolar type 2",
                       "17"="Alveolar type 2",
                       "18"="Mast cells",
                       "19"="Fibroblasts:CCL19",
                       "21"="Smooth muscle:Airway",
                       "22"="Hybrid cells",
                       "24"="Endothelial cells:Smooth muscle",
                       "25"="Alveolar type 2",
                       "26"="Ciliated cells:2",
                       "27"="Endothelial cells:Lymphatic",
                       "28"="Pre-ciliated cells",
                       "29"="Secretory cells",
                       "30"="B cells:Plasma",
                       "31"="Unknown",
                       "32"="Ciliated cells:3",
                       "33"="Ionocytes/NEC",
                       "34"="Alveolar type 1",
                       "35"="Neurons",
                       "36"="Secretory cells",
                       "37"="Cartilage",
                       "38"="Unknown"
)
VU19[["labels"]] = as.character(Idents(VU19))
meta.data2[rownames(VU19@meta.data),"labels"]  = VU19$labels
table(meta.data2$labels)
meta.data2$labels %<>% gsub("Ciliated intermediate cells","Unknown",.)

object@meta.data = meta.data2
Idents(object) <- "labels"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "labels", colors = Singler.colors)
table(object$conditions)
UMAPPlot.1(object, cols = ExtractMetaColor(object),
           label = F, label.repel = T,pt.size = 0.3,
           label.size = 4, repel = T,no.legend = T,do.print = T,
           do.return = F,alpha = 0.9,
           title = "Cell types",unique.name = "conditions")


save(object,file="data/Lung_16_distal_20191022.Rda")

DefaultAssay(object) <- 'RNA'
for(i in c(5,9,10,12,15)/10){
        object %<>% FindClusters(resolution = i)
        Idents(object) = paste0("RNA_snn_res.",i)
}
# follow Renat's instruction : "16 w COPD VU-19 clusters"
object@meta.data$labels %<>% as.character()
object@meta.data$labels = object@meta.data$labels.copy
object@meta.data$labels.copy <-  as.character(object@meta.data$labels)
object@meta.data[object$RNA_snn_res.0.5 == 8,"labels"] = "Alveolar type 2"
object@meta.data[object$RNA_snn_res.0.5 == 24,"labels"] = "Alveolar type 1"
object@meta.data[object$RNA_snn_res.0.5 == 12 & 
                         object$labels != "Dendritic cells:Plasmacytoid","labels"] = "B cells"
object@meta.data[object$RNA_snn_res.0.5 == 14,"labels"] = "Mast cells"
object@meta.data[object$RNA_snn_res.0.5 == 15,"labels"] = "Neutrophils"
object@meta.data[object$RNA_snn_res.0.5 == 16,"labels"] = "Hybrid cells"
object@meta.data[object$RNA_snn_res.0.5 == 18,"labels"] = "Endothelial cells:Lymphatic"
object@meta.data[object$RNA_snn_res.0.5 == 20,"labels"] = "Pre-ciliated cells"
object@meta.data[object$RNA_snn_res.0.5 == 22,"labels"] = "B cells:Plasma"
object@meta.data[object$RNA_snn_res.0.5 == 25,"labels"] = "Neurons"
object@meta.data[object$RNA_snn_res.0.5 == 27,"labels"] = "Cartilage"
object@meta.data[object$RNA_snn_res.0.5 == 26,"labels"] = "Serous cells"
object@meta.data[object$RNA_snn_res.0.5 == 28,"labels"] = "Red blood cells"

object@meta.data[object$RNA_snn_res.0.8 == 1 & object$orig.ident == "VU19-D" ,
                 "labels"] = "Ciliated cells:1"
object@meta.data[object$RNA_snn_res.0.8 == 9 & object$orig.ident == "VU19-D" ,
                 "labels"] = "Ciliated cells:2"
object@meta.data[object$RNA_snn_res.0.8 == 14 & object$orig.ident == "VU19-D" ,
                 "labels"] = "Ciliated cells:3"
object@meta.data[object$RNA_snn_res.0.8 == 26 & object$orig.ident == "VU19-D" ,
                 "labels"] = "Ciliated cells:2"
object@meta.data[object$RNA_snn_res.0.8 == 32,"labels"] = "Hybrid cells"

object@meta.data[object$RNA_snn_res.0.9 == 3,"labels"] = "Endothelial cells:HEVs"
object@meta.data[object$RNA_snn_res.0.9 == 5,"labels"] = "Endothelial cells:Capillary"
object@meta.data[object$RNA_snn_res.0.9 == 25,"labels"] = "Endothelial cells:Smooth muscle"
object@meta.data[object$RNA_snn_res.0.9 == 8,"labels"] = "Secretory cells"
object@meta.data[object$RNA_snn_res.0.9 == 9,"labels"] = "Basal cells"
object@meta.data[object$RNA_snn_res.0.9 %in% c(24,30),"labels"] = "Intermediate cells"
#object@meta.data[object$RNA_snn_res.1 %in% c(1,14,19,20,34),"labels"] = "Ciliated cells:1"
#object@meta.data[object$RNA_snn_res.1 %in% 23,"labels"] = "Ciliated cells:1"

object@meta.data[object$RNA_snn_res.1.2 == 22,"labels"] = "Fibroblasts:IGF1"
object@meta.data[object$RNA_snn_res.1.2 == 4,"labels"] = "Fibroblasts:TBX4"
object@meta.data[object$RNA_snn_res.1.2 == 18,"labels"] = "Fibroblasts:CCL19"

object@meta.data[object$RNA_snn_res.1.5 == 5,"labels"] = "Smooth muscle:Vascular"
#object@meta.data[object$RNA_snn_res.1.5 == 26,"labels"] = "Smooth muscle:Airway"
object@meta.data[object$RNA_snn_res.1.5 == 39,"labels"] = "Pericytes"
object@meta.data[object$RNA_snn_res.1.5 %in% c(1,3),"labels"] = "T cells:TRM"
object@meta.data[object$RNA_snn_res.1.5 == 7,"labels"] = "T cells:CD4+"
object@meta.data[object$RNA_snn_res.1.5 == 10,"labels"] = "NK cells"
#Idents(object) = "RNA_snn_res.1.5"
#cluster_43 <- subset(object, idents = 43)
#FeaturePlot.1(cluster_43, features = c("CD163","MKI67"),border = T)
#Mf_proliferating  <- subset(cluster_43, subset = `CD163` > 0 & `MKI67` >0)
object@meta.data[colnames(Mf_proliferating),"labels"] = "Macrophages:Proliferating"
cluster_43@meta.data = cbind(cluster_43@meta.data, cluster_43@reductions$umap@cell.embeddings)

object@meta.data[colnames(cluster_43)[cluster_43$UMAP_1 < -6],"labels"] = "T cells:Proliferating"
object@meta.data[object@assays$RNA@data["GJA5",]>0 & 
                         object$RNA_snn_res.0.9 %in% c(3,5),"labels"] = "Endothelial cells:Arterial"
object@meta.data[object@assays$RNA@data["SCGB3A2",]>1 & 
                         object$RNA_snn_res.0.9 %in% c(8,24,30),"labels"] = "Secretory cells:Distal"
object@meta.data[object@assays$RNA@data["TFF3",]>4 & 
                         object$RNA_snn_res.0.9 %in% c(8,24),"labels"] = "Mucus-producing:Goblet"
object@meta.data$major_labels = gsub(":.*","",as.character(object$labels))

save(object, file = "data/Lung_16_distal_20191022.Rda")



table(object$conditions)
COPD_cells <- colnames(object)[object$conditions %in% "COPD"]
distal_cells <- colnames(object)[object$conditions %in% "distal"]
distal_cells <- distal_cells[sample(1:46400, 15145)]
subset_object <- object[,c(COPD_cells,distal_cells)]

UMAPPlot.1(subset_object, group.by = "conditions",
           label = F, label.repel = T,pt.size = 0.3,
           label.size = 2, repel = T,no.legend = F,do.print = T,
           do.return = F,alpha = 0.9,
           title = "distal vs COPD",unique.name = "conditions")

