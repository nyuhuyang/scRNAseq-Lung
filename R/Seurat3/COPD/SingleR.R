library(SingleR)
library(Seurat)
library(scMerge)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file="data/Lung_24distal_20191004.Rda"))
DefaultAssay(object) = "SCT"
object_distal <- object

(load(file="data/Lung_24terminal_20191004.Rda"))
DefaultAssay(object) = "SCT"
object_terminal <- object

(load(file="data/Lung_24proximal_20191004.Rda"))
DefaultAssay(object) = "SCT"
object_proximal <- object

object <- Reduce(function(x, y) merge(x, y, do.normalize = F), 
                 list(object_distal,object_terminal,object_proximal))
table(object$cell.types)

remove(object_distal,object_terminal,object_proximal);GC()

sec_lung <- as.SingleCellExperiment(object)

(load(file="data/Lung_16_distal_20191017.Rda"))
DefaultAssay(object) = "SCT"
sce <- as.SingleCellExperiment(object)

singler = SingleR(test = sce, ref = sec_lung, 
                  labels = sec_lung$cell.types)

save(singler,file="output/singler_Lung_16_distal_20191017.Rda")
