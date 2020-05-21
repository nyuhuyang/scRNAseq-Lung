########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load files
object = readRDS(file = "data/Lung_28_Global_20200511.rds") 
# Need 64GB
DefaultAssay(object) = "SCT"
Idents(object) = "annotations"
cell.types <- unique(Idents(object))
if(any(grep("unknown",cell.types))) cell.types = cell.types[-grep("unknown",cell.types)]
cell.type = as.character(cell.types[args])
system.time(Lung_markers <- FindMarkers.UMI(object, ident.1 = cell.type, ident.2 = NULL,
                                               logfc.threshold = 0.1, only.pos = T,
                                               test.use = "MAST",min.cells.group = 2))
Lung_markers$gene = rownames(Lung_markers)
Lung_markers$cluster = cell.type
if(args < 10) args = paste0("0", args)
write.csv(Lung_markers,paste0(path,"Lung_28-",args,"_FC0.1_",cell.type,".csv"))