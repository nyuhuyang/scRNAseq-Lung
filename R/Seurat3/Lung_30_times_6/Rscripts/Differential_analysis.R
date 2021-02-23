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
plan("multiprocess", workers = 4)
plan()
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

step = 3 # DEGs for cell types at d28 only
if(step == 1){ # need 32 GB
        opts = data.frame(resolution = c(rep(0.4,16),
                                         rep(1,21),
                                         rep(2,34)),
                          cluster = c(0:15,
                                      0:20,
                                      0:33),
                          stringsAsFactors = F)
        set.seed(101)
        print(res <- opts$resolution[args])
        print(cluster <- opts$cluster[args])
        
        object = readRDS(file = "data/Lung_time_6_20201001.rds")
        DefaultAssay(object)  = "SCT"
        object %<>% FindClusters(resolution = res)
        Idents(object) = paste0("SCT_snn_res.",res)
        system.time(Lung_markers <- FindMarkers.UMI(object, 
                                                    ident.1 = cluster,
                                                    latent.vars = "nFeature_SCT",
                                                    logfc.threshold = 0.05,
                                                    only.pos = T,
                                                    test.use = "MAST"))
        Lung_markers$cluster = cluster
        Lung_markers$gene = rownames(Lung_markers)
        write.csv(Lung_markers,paste0(path,"Lung_time_6_FC1-res=",res,"_cluster=",cluster,".csv"))
}

if(step == 2){ # need 64 GB # find marker genes for each cell types
        object = readRDS(file = "data/Lung_time_6_20201001.rds")
        
        Idents(object) = "Doublets"
        object %<>% subset(idents = "Singlet")
        Idents(object) = "annotations3"
        clusters =  c("BC-0","BC-m","BC-p","BC-S","BC1","BC2",
                      "C1","C2","e-S","IC-S","IC1","IC2","IC2-s",
                      "Ion","NEC","p-C","S","S-d","Un")

        print(cluster <- clusters[args])
        
        DefaultAssay(object)  = "SCT"

        system.time(Lung_markers <- FindMarkers.UMI(object, 
                                                    ident.1 = cluster,
                                                    latent.vars = "nFeature_SCT",
                                                    logfc.threshold = 0.05,
                                                    only.pos = F,
                                                    test.use = "MAST"))
        Lung_markers$cluster = cluster
        Lung_markers$gene = rownames(Lung_markers)
        write.csv(Lung_markers,paste0(path,"Lung_time_6_FC0.05_annotation=",cluster,".csv"))
}

if(step == 3){ # need 64 GB # DEGs for cell types at d28 only
        object = readRDS(file = "data/Lung_time_6_20201001.rds")
        
        Idents(object) = "Doublets"
        object %<>% subset(idents = "Singlet")
        Idents(object) = "conditions"
        object %<>% subset(idents = "Day-28")
        
        Idents(object) = "annotations3"
        clusters =  c("BC-m","BC-p","BC2","C1","C2","e-S","IC-S","IC1", "IC2",
                      "IC2-s","Ion","NEC","p-C", "S","S-d")
        
        print(cluster <- clusters[args])
        
        DefaultAssay(object)  = "SCT"
        
        system.time(Lung_markers <- FindMarkers.UMI(object, 
                                                    ident.1 = cluster,
                                                    latent.vars = "nFeature_SCT",
                                                    logfc.threshold = 0.05,
                                                    only.pos = F,
                                                    test.use = "MAST"))
        Lung_markers$cluster = cluster
        Lung_markers$gene = rownames(Lung_markers)
        write.csv(Lung_markers,paste0(path,"Lung_time_6_FC0.05_Day-28_annotation=",cluster,".csv"))
}
