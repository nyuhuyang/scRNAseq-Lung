########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# Need 64GB
# load files
object = readRDS(file = "data/Lung_30_20200710.rds") 
# Need 64GB
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
Idents(object) = "annotations3"
object <- subset(object, idents = "Un", invert = TRUE)
cell.types <- sort(as.character(unique(Idents(object))))
cell.type = cell.types[args]
system.time(Lung_markers <- FindMarkers.UMI(object, ident.1 = cell.type, ident.2 = NULL,
                                               logfc.threshold = 0.1, only.pos = T,
                                               test.use = "MAST",min.cells.group = 2))
Lung_markers$gene = rownames(Lung_markers)
Lung_markers$cluster = cell.type
if(args < 10) args = paste0("0", args)
write.csv(Lung_markers,paste0(path,"Lung_30-",args,"_FC0.1_",cell.type,".csv"))