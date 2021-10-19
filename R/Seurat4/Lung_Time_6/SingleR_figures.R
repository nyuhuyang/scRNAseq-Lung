# conda activate r4.0.3
library(Seurat)
library(magrittr)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
object = readRDS(file = "data/Lung_time6_20210908.rds")
pred = readRDS(file = "output/Lung_time6_20210908_Lung30_pred.rds")

singlerDF = data.frame("Cell_subtype" = pred$pruned.labels,
                       row.names = rownames(pred))
table(rownames(pred) == rownames(object@meta.data))
table(is.na(singlerDF$Cell_subtype))
singlerDF$Cell_subtype[is.na(singlerDF$Cell_subtype)]= "Un"

##############################
# adjust cell label
##############################
# combine cell types

##############################
# process color scheme
##############################
table(colnames(object) == rownames(singlerDF))
object$Cell_subtype.lung30 = singlerDF[,"Cell_subtype"]

meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")

meta.data = meta.data[!duplicated(meta.data$Cell_subtype),]
meta.data = meta.data[order(meta.data$Cell_subtype),]
object$Cell_subtype.lung30.colors = plyr::mapvalues(object$Cell_subtype.lung30,
                                                    from = meta.data$Cell_subtype,
                                                    to = meta.data$Cell_subtype.colors)
    
UMAPPlot.1(object = object, label = T, label.repel = T,group.by = "Cell_subtype",
    no.legend = T,#cols = Singler.colors,
    pt.size = 0.1,label.size = 5,alpha = 0.85,
    do.print = T,do.return = F,
    title ="labeling by Lung 30")

saveRDS(meta.data, "output/20211016/meta.data_Cell_subtype_time6.rds")

