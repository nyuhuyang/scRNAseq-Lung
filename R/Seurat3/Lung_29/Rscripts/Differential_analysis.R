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

step = 2
if(step == 1){ # need 128 GB
        opts = data.frame(cluster = c(0:76),
                          stringsAsFactors = F)
        set.seed(101)
        (cluster = opts$cluster[args])
        
        object = readRDS(file = "data/Lung_29_20200617.rds")
        DefaultAssay(object)  = "SCT"
        object %<>% FindClusters(resolution = 2)
        Idents(object) = "SCT_snn_res.2"
        system.time(Lung_markers <- FindMarkers.UMI(object, 
                                                    ident.1 = cluster,
                                                    latent.vars = "nFeature_SCT",
                                                    logfc.threshold = 0, 
                                                    only.pos = T,
                                                    test.use = "MAST"))
        Lung_markers$cluster = cluster
        Lung_markers$gene = rownames(Lung_markers)
        write.csv(Lung_markers,paste0(path,"Lung_29-res=2_cluster=",cluster,".csv"))
}

if(step == 2){ # need 128 GB
        object = readRDS(file = "data/Lung_29_20200617.rds")
        Idents(object) = "annotations3"
        cell.types <- unique(Idents(object))
        cell.type = as.character(cell.types[args])
        
        system.time(Lung_markers <- FindAllMarkers.UMI(object, 
                                                       ident.1 = cell.type,
                                                       ident.2 = NULL,
                                                       latent.vars = "nFeature_SCT",
                                                       logfc.threshold = 0.1, only.pos = T,
                                                       test.use = "MAST",min.cells.group = 2))
        Lung_markers$gene = rownames(Lung_markers)
        Lung_markers$cluster = cell.type
        if(args < 10) args = paste0("0", args)
        write.csv(Lung_markers,paste0(path,"Lung_29-",args,"_FC0.1_",cell.type,".csv"))
}