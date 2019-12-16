########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
library(MAST)
library(future)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# change the current plan to access parallelization
plan("multiprocess", workers = 8)
plan()

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# load data
(load(file = paste0("data/Lung_24_20191206.Rda")))
DefaultAssay(object)  = "SCT"
Idents(object) = "cell.types"
object %<>% sortIdent
cell.types <- unique(Idents(object))
cell.type = cell.types[args]
print(paste("FindMarkers for=",cell.type))
sub_object <- subset(object, idents = cell.type)
Lung_markers <- FindAllMarkers.UMI(sub_object, logfc.threshold = 0, only.pos = F,
                                   return.thresh = 1,
                                   test.use = "MAST")
write.csv(Lung_markers,paste0(path,"Lung_24-FC0_markers_",args,"_",cell.type,".csv"))
# Differential analysis
#Top_n = 5
#top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
#object %<>% ScaleData(features=unique(c(as.character(top$gene))))

#DoHeatmap.1(object, marker_df = Lung_markers, Top_n = 5, do.print=T,
#            angle = 0,group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
#            assay = "SCT",
#            label=T, cex.row=5, legend.size = 5,width=10, height=7,unique.name = T,
#            title = paste("Top 5 markers of each clusters in 24 sampels"))

#res = read.csv("~/Downloads/Lung_24-FC0.75_markers_res=0.8.csv")
#res %>% group_by(cluster) %>% top_n(20,avg_logFC) %>% .["avg_logFC"] %>% min()
