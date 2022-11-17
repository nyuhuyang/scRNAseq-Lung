#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
invisible(lapply(c("Seurat","SingleR","SingleCellExperiment","magrittr","data.table","Matrix"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

test_df = data.frame("datasets" = rep(c("WholeLung","HighResolution"),each =2),
                     "references" = rep(c("scRNA-Lung30","bulkRNA-Lung30"), time = 2))

# ====== load single cell =============
dataset <- test_df$datasets[args]
object <- readRDS(paste0("data/GSE141259_",dataset,".rds"))
sce <- SingleCellExperiment(list(logcounts=object[["RNA"]]@data),
                            colData=DataFrame(object@meta.data))
rm(object);GC()

# ====== load Lung30 data =============
reference <- test_df$references[args]
Lung = readRDS(file = "data/Lung_SCT_30_20210831.rds")
if(reference == "scRNA-Lung30"){
    sce_Lung <- SingleCellExperiment(list(logcounts=Lung[["SCT"]]@data),
                                colData=DataFrame(Lung@meta.data))
} else if(reference == "bulkRNA-Lung30"){
    Lung$Cell_subtype_orig_ident <- paste0(Lung$Cell_subtype,"_",Lung$orig.ident)
    Cell_subtype_orig_ident_exp <- AverageExpression(Lung,assays = "SCT",group.by = "Cell_subtype_orig_ident")
    sce_Lung <- SingleCellExperiment(list(logcounts=Cell_subtype_orig_ident_exp$SCT),
                                                   colData=DataFrame("Cell_subtype" = gsub("_.*","",colnames(Cell_subtype_orig_ident_exp$SCT)),
                                                                     row.names = colnames(Cell_subtype_orig_ident_exp$SCT)))
    }
rownames(sce) %<>% toupper()
common <- Reduce(intersect, list(rownames(sce),
                                 rownames(sce_Lung)
))
length(common)
table(sce_Lung$Cell_subtype)
system.time(trained <- trainSingleR(ref = sce_Lung[common,],
                                    labels=sce_Lung$Cell_subtype))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = paste0("output/GSE141259_",dataset,"_singelR_",reference,".rds"))
