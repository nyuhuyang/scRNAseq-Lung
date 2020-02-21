########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
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

cell.types = c("AT-p","AT1","AT2","B","BC","BC-p",
               "C1","C2","C3","C4","Cr","DC",
               "DC-p","En-A","En-C","En-L","En-SM","En-V",
               "F1","F2","F3","F4","F5","Gli",
               "H","IC","Ion","M-p","M0","M1",
               "M2","MC","MEC","Mon","NEC","Neu",
               "Nr","p-C","PC","Pr","S","S-d",
               "SM1","SM2","SM3","SMG-Muc","SMG-Ser","Sq",
               "T-7s","T-cn","T-inf","T-NK","T-p","T-reg","T-rm") 
(c = cell.types[args])
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load files
object = readRDS(file = "data/Lung_28_Global_20200206.rds") 

Idents(object) = "cell.labels"
system.time(Lung_markers <- FindMarkers.UMI(object, ident.1 = c, ident.2 = NULL,
                                               logfc.threshold = 0.5, only.pos = T,
                                               test.use = "MAST",min.cells.group = 2))
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
if(args < 10) args = paste0("0", args)
write.csv(Lung_markers,paste0(path,"Lung_28-",args,"_",c,".csv"))
