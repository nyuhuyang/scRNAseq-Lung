library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#====== 2.1 Identify cell types ==========================================
(load(file="data/Lung_28_20200116.Rda"))
#======== rename ident =================
Idents(object) = "SCT_snn_res.0.8"
object %<>% RenameIdents("0" = "T",
    "1" = "En",
    "2" = "F",
    "3" = "T",
    "4" = "En",
    "5" = "Mon",
    "6" = "SAE-H",
    "7" = "T",
    "8" = "Neu",
    "9" = "B",
    "10" = "SM",
    "11" = "SAE-Sec1",
    "12" = "SAE-C1",
    "13" = "SAE",
    "14" = "SAE-BC",
    "15" = "AT2",
    "16" = "Mac",
    "17" = "DC",
    "18" = "D",
    "19" = "AT2",
    "20" = "Mac",
    "21" = "SAE-Ion",
    "22" = "MC",
    "23" = "SAE-IC1",
    "24" = "SMG",
    "25" = "SAE-IC2",
    "26" = "SM",
    "27" = "T",
    "28" = "En-Lym",
    "29" = "SM3",
    "30" = "Per",
    "31" = "P",
    "32" = "SAE-C",
    "33" = "SAE-C-pre",
    "34" = "Car",
    "35" = "PC",
    "36" = "SAE-C",
    "37" = "SAE-Sq",
    "38" = "T",
    "39" = "NEC",
    "40" = "AT1",
    "41" = "Nr",
    "42" = "SAE",
    "43" = "SAE")
object[["cell_types"]] = as.character(Idents(object))
meta.data <- cbind(object@meta.data,object@reductions$umap@cell.embeddings)
object@meta.data[(meta.data$UMAP_1 > 5 & meta.data$UMAP_2 < -5),
                 "cell_types"] = "AT2"
object$cell_types %<>% gsub("-.*","",.) %>% gsub("[0-9]+","",.)

Idents(object) = "cell_types"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "cell_types", colors = c(Singler.colors,Singler.colors))

TSNEPlot.1(object, group.by = "cell_types",cols = ExtractMetaColor(object),label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "Cell types in all 28 samples")
UMAPPlot.1(object, group.by = "cell_types",cols = ExtractMetaColor(object),label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "Cell types in all 28 samples")


conditions = c("proximal", "distal","terminal","COPD")
Idents(object) = "conditions"
for(i in seq_along(conditions)){
    sub_object <- subset(object, idents = conditions[i])
    Idents(sub_object) = "SCINA"
    
    UMAPPlot.1(sub_object, group.by = "SCINA",cols = ExtractMetaColor(sub_object),label = T,
               label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
               unique.name = "conditions",
               do.print = T,do.return = F,
               title = paste("Cell types in",conditions[i]))
    TSNEPlot.1(sub_object, group.by = "cell_types",cols = ExtractMetaColor(sub_object),label = T,
               label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
               unique.name = "conditions",
               do.print = T,do.return = F,
               title = paste("Cell types in",conditions[i]))
    Progress(i,length(conditions))
}

write.csv(table(object$cell_types,object$SCINA),
          paste0(path,"Inherited_vs_SCINA.csv"))