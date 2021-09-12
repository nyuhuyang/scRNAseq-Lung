library(Seurat)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

opts = paste0("SCT_snn_res.",c(0.01, 0.1, 0.2, 0.5))
opt = opts[args]
print(opt)
#==========================
object = readRDS(file = "data/Lung_RNA_62_20210831.rds")

Idents(object) = opt

markers = FindAllMarkers_UMI(object,
                          group.by = opt,
                          assay = "RNA",
                          logfc.threshold = 0.5,
                             only.pos = T,
                             test.use = "MAST",
                             latent.vars = "nFeature_SCT")

write.csv(markers,paste0(path,args,"_",opt, ".csv"))
