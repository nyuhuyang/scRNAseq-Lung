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
 object %<>% FindClusters(resolution = 0.06)
Idents(object) = "SCT_snn_res.0.06"
object %<>% RenameIdents("0" = "T",
                         "1" = "Myeloid",
                         "2" = "En",
                         "3" = "SAE",
                         "4" = "SAE",
                         "5" = "F",
                         "6" = "AT",
                         "7" = "SM",
                         "8" = "B",
                         "9" = "MC",
                         "10" = "SMG",
                         "11" = "En-Lym",
                         "12" = "NEC_Car_Nr")
object[["cell.types"]] = as.character(Idents(object))
Idents(object) = "cell.types"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "cell.types", colors = Singler.colors)

lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun) 
        fun(object, group.by = "cell.types",cols = ExtractMetaColor(object),
            label = T,
            label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
            do.print = T,do.return = F,title = "Cell types in all 28 samples"))
save(object, file = "data/Lung_28_harmony_20200131.Rda")
