########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","tidyr","magrittr","gplots","MAST",
                   "future"), function(x) {
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
# require 32GB
# load data
object = readRDS(file = "data/Lung_28_Global_20200511.rds")
DefaultAssay(object)  = "SCT"
Idents(object) = "annotations"
cell.types <- unique(Idents(object))
cell.types = cell.types[-grep("unknown",cell.types)]
cell.type = cell.types[args]
print(paste("FindMarkers for=",cell.type))
Lung_markers <- FindMarkers.UMI(object, ident.1 = cell.type, 
                                logfc.threshold = 0, only.pos = F,
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
