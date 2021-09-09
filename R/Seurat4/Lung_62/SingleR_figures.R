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
pred = readRDS("output/Lung_SCT_59_20210814_singleR_pred.rds")
object = readRDS("data/Lung_SCT_59_20210814.rds")


singlerDF = data.frame("cell_types" = gsub("\\.","-",pred$pruned.labels),
                       row.names = rownames(pred))

table(is.na(singlerDF$cell_types))
singlerDF$cell_types[is.na(singlerDF$cell_types)]= "Un"


##############################
# adjust cell color
##############################

# combine cell types
Lung = readRDS(file = "data/Lung_SCT_30_20200710.rds")
meta.data = Lung@meta.data
#rm(Lung);GC()
meta.data = meta.data[!duplicated(meta.data$cell_types),]
meta.data = meta.data[order(meta.data$cell_types),]
duplicated_color = meta.data$cell_types.colors[duplicated(meta.data$cell_types.colors)]
meta.data$cell_types.colors[duplicated(meta.data$cell_types.colors)] = c("#ffb399","#ffd34d","#c06000","#08b000","#0078a4","#e6ad00","#ff86bb")
singlerDF$cell_types.colors = singlerDF$cell_types
singlerDF$cell_types.colors %<>% plyr::mapvalues(from = meta.data$cell_types,
                                                 to = meta.data$cell_types.colors)
singlerDF$cell_types.colors %<>% plyr::mapvalues(from = c("BC","SM1-2"),
                                                 to = c("#70AD47","#C8D8E0"))
object@meta.data %<>% cbind(singlerDF)

lapply(c(TSNEPlot.1,UMAPPlot.1), function(fun)
    fun(object = object, label = T, label.repel = T,group.by = "cell_types",
        no.legend = T,
        pt.size = 0.1,label.size = 3,
        do.print = T,do.return = F,
        title ="labeling by Lung 30 scRNA-seq"))

saveRDS(object, file = "data/Lung_SCT_59_20210814.rds")
