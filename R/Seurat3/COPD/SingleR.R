library(SingleR)
library(Seurat)
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

sec_distal <- as.SingleCellExperiment(object)
system.time(trained <- trainSingleR(sec_distal,labels = sec_distal$cell.types))

(load(file="data/Lung_16_distal_20191017.Rda"))
DefaultAssay(object) = "SCT"
sce <- as.SingleCellExperiment(object)
remove(object);GC()

singler = SingleR(test = sce, ref = sec_lung,  fine.tune = FALSE,
                  labels = sec_lung$cell.types)

save(singler,file="output/singlerF_Lung_16_distal_20191017.Rda")
