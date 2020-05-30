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
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# change the current plan to access parallelization

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))
opts = data.frame(resolution = c(rep(1,47),
                                 rep(2,64),
                                 rep(3,84),
                                 rep(4,102),
                                 rep(5,119)),
                  cluster = c(0:46,
                              0:63,
                              0:83,
                              0:101,
                              0:118),
                  stringsAsFactors = F)
set.seed(101)
(res = opts$resolution[args])
(cluster = opts$cluster[args])
step = 1
if(step == 1){ # need 128 GB
        (load(file = "data/Lung_28_Nointeg_20200131.Rda"))
        DefaultAssay(object)  = "SCT"
        object %<>% FindClusters(resolution = res)
        Idents(object) = paste0("SCT_snn_res.",res)
        system.time(Lung_markers <- FindMarkers.UMI(object, 
                                                    ident.1 = cluster,
                                                    latent.vars = "nFeature_SCT",
                                                    logfc.threshold = 0, 
                                                    only.pos = T,
                                                    test.use = "MAST"))
        Lung_markers$cluster = cluster
        Lung_markers$gene = rownames(Lung_markers)
        write.csv(Lung_markers,paste0(path,"Lung_28-res=",res,"_cluster=",cluster,".csv"))
}
