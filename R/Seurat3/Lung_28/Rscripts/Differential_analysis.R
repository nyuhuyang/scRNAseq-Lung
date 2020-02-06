########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr",
                   "MAST","future","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# change the current plan to access parallelization

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))
res = list(c(1.9,1),
           c(1.9,0.5),
           c(1.9,0.25),
           c(2.7,1),
           c(2.7,0.5),
           c(2.7,0.25))
(r <- res[[args]])
step = 1
if(step == 1){ # need 128 GB
        (load(file = "data/Lung_28_harmony_20200131.Rda"))
        DefaultAssay(object)  = "SCT"
        object %<>% FindClusters(resolution = r[1])
        Idents(object) = paste0("SCT_snn_res.",r[1])
        system.time(Lung_markers <- FindAllMarkers.UMI(object, 
                                                       logfc.threshold = r[2], only.pos = T,
                                           test.use = "MAST"))
        Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
        write.csv(Lung_markers,paste0(path,"Lung_28-FC",r[2],"_res=",r[1],".csv"))
}
