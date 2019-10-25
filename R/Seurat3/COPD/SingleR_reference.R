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

Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
anyNA(object$cell.types)
Idents(object) = "cell.types"
object %<>% sortIdent()
object <- subset(object, idents =c("Basal-like cells",
                                   "Basal cells:Proliferating",
                                   "Ciliated intermediate cells",
                                   "Macrophages-like",
                                   "Mucous gland cells",
                                   "Myo-epithelial cells",
                                   "Proliferating basal cells",
                                   "Serous gland cells",
                                   "T cells:7SK.2+",
                                   "T cells:7SK.2",
                                   "Unknown"), invert = T)
object@meta.data$cell.types %<>% as.character()
object$cell.types %<>% gsub("^Distal secretory cells$","Secretory cells:Distal",.)
object$cell.types %<>% gsub("^Endothelial cells:arterial$","Endothelial cells:Arterial",.)
object$cell.types %<>% gsub("^Endothelial cells:HEV$","Endothelial cells:HEVs",.)
object$cell.types %<>% gsub("^Endothelial cells:smooth muscel$","Endothelial cells:Smooth muscle",.)
object$cell.types %<>% gsub("^Endothelial:Smooth muscle$","Endothelial cells:Smooth muscle",.)
object$cell.types %<>% gsub("^Goblet$","Secretory cells:Goblet",.)
object$cell.types %<>% gsub("^Hybrid$","Hybrid cells",.)
object$cell.types %<>% gsub("^Macrophage.*","Macrophages",.)
object$cell.types %<>% gsub("^Neuro endocrine$","Neuroendocrine",.)
object$cell.types %<>% gsub("^Neutrophil$","Neutrophils",.)
object$cell.types %<>% gsub("^Plasma cells$","B cells:Plasma",.)

sec_lung <- as.SingleCellExperiment(object)
#save(sec_lung,file="output/sec_lung_24_20191024.Rda")
#(load(file="output/sec_lung_24_20191024.Rda"))
object <- as.Seurat(sec_lung)
Idents(object) = "cell.types"
object %<>% sortIdent()
Lung_exp <- AverageExpression(object = object, assays = "RNA")

# for new SingleR ======================
Lung_exp$RNA
Lung_sce <- SingleCellExperiment(
        assays = list(logcounts = Lung_exp$RNA), 
        colData = list(labels = colnames(Lung_exp$RNA))
)

system.time(trained <- trainSingleR(Lung_exp$RNA,labels = colnames(Lung_exp$RNA)))

save(trained,file="output/singlerT_Lung_16_distal_20191023.Rda")

# for old SingleR ======================
# Create Singler Reference
ref = CreateVariableGeneSet(ref_data = Lung_exp$RNA,
                            types = colnames(Lung_exp$RNA),
                            n = 500)
save(ref,file="data/ref_Lung_24_20191024.Rda")
