########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","magrittr","MAST",
                   "future","ggplot2","tidyr","harmony"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))
step = 1
if(step == 1){ # 32GB
        df_samples <- readxl::read_excel("doc/Annotations-4-30-20-RS.xlsx",
                                 sheet = "Sheet2")
        (sample = df_samples[args,])
        object = readRDS(file = paste0("output/Lung_28_",sample$`Group #`,"_",
                                      sample$`Group name`,"_",
                                      sample$`Method`,"-2020406.rds"))
        DefaultAssay(object) = "SCT"
        object %<>% FindClusters(resolution = sample$Resolution)
        object %<>% subset(idents = as.integer(strsplit(sample$Cluster,",")[[1]]))
        saveRDS(colnames(object), file = paste0(path, "Annotation-",sample$`Cell type`,".rds"))
}