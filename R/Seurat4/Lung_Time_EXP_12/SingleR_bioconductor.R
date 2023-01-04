#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)
library(data.table)
library(Matrix)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# ======= load single cell dataset =================
object = readRDS(file = "data/Lung_time6_20210908.rds")

sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ====== load lung 30 data =============
Lung = readRDS(file = "data/Lung_SCT_30_20210831.rds")
sce_lung <- SingleCellExperiment(list(logcounts=Lung[["SCT"]]@data),
                            colData=DataFrame(Lung@meta.data))
rm(Lung);GC()

#============
common <- Reduce(intersect, list(rownames(sce),
                                 rownames(sce_lung)
))
length(common)
table(sce_lung$Cell_subtype)
system.time(trained <- trainSingleR(ref = sce_lung[common,],
                                    labels=sce_lung$Cell_subtype))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/Lung_time6_20210908_Lung30_pred.rds")
