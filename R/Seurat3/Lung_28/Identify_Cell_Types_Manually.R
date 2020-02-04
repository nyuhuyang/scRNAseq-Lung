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
(load(file="data/Lung_28_harmony_20200131.Rda"))
#======== rename ident =================
(res = c(0.004,0.01))
for(i in seq_along(res)) object %<>% FindClusters(resolution = res[i])
Idents(object) = "SCT_snn_res.0.004"
object %<>% RenameIdents("0" = "T",
                         "1" = "En_&_SM",
                         "2" = "Myeloid",
                         "3" = "Epi",
                         "4" = "Epi",
                         "5" = "F_&_Nr",
                         "6" = "Epi",
                         "7" = "B",
                         "8" = "MC",
                         "9" = "En-Lym")
object[["cell.types"]] = as.character(Idents(object))
Idents(object) = "cell.types"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "cell.types", colors = Singler.colors)

TSNEPlot.1(object, group.by = "cell.types",cols = ExtractMetaColor(object),
           label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "Cell types in all 28 samples")
UMAPPlot.1(object, group.by = "cell.types",cols = ExtractMetaColor(object),
           label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "Cell types in all 28 samples")
save(object, file = "data/Lung_28_harmony_20200131.Rda")
