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
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
set.seed=101
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))


# samples
samples = c("All","Day-0","Day-3","Day-7","Day-14","Day-21","Day-28",
             "Day-56","Day-122")
(sample = samples[args])

# 3.1.1 load data
(load(file = "data/Lung_8_time_20190808.Rda"))
Idents(object) <-  "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) <- "orig.ident"
if(args>1) object %<>% subset(idents = sample)
Idents(object) <- "integrated_snn_res.0.6"
Lung_markers <- FindAllMarkers.UMI(object, logfc.threshold = 0.1,
                                   only.pos = T)

write.csv(Lung_markers,paste0(path,"Lung_8-",sample,"_markers.csv"))
#Lung.markers.csv =read.csv(file = paste0(path,"Lung.markers.csv")

#object %<>% ScaleData(features=unique(Lung_markers$gene))
#DoHeatmap.1(object, marker_df = Lung_markers, Top_n = 5, do.print=T, angle = 0,
#            group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
#            assay = "SCT",
#            label=T, cex.row=5, legend.size = 5,width=10, height=7,unique.name = T,
#            title = paste("Top 5 markers in each clusters",sample,
#                          ifelse(sample=="All","sampels","sample")))

