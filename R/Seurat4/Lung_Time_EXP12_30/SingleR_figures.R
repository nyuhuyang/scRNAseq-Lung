# conda activate r4.1.3
library(Seurat)
library(magrittr)
#library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
library(data.table)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame for Lung 30
###############################
pred1 <- readRDS("output/Lung_time15_20220523_pseduBulk_human_lung_v2_singleR_pred.rds")
pred2 <- readRDS("output/Lung_time15_20220523_pseduBulk_pseduBulk_Lung30_singleR_pred.rds")

meta.data <- readRDS("output/Lung_time15_metadata_20220523_v2.rds")

table(rownames(pred1) == rownames(pred2))
singlerDF = data.frame("ann_finest_level" = pred1$pruned.labels,
                       "Cell_subtype" = pred2$pruned.labels,
                       row.names = rownames(pred1))

table(is.na(singlerDF$ann_finest_level))
table(is.na(singlerDF$Cell_subtype))

singlerDF$ann_finest_level[is.na(singlerDF$ann_finest_level)]= "unknown"
singlerDF$Cell_subtype[is.na(singlerDF$Cell_subtype)]= "unknown"

table(singlerDF$ann_finest_level)
table(singlerDF$Cell_subtype)
##############################
# adjust cell label
##############################
# combine cell types
metadata <- readRDS("output/Lung_30_20210831_metadata_v2.rds")
metadata = metadata[!duplicated(metadata$Cell_subtype),c("Cell_subtype","Cell_subtype.colors")]

singlerDF$Cell_subtype.colors = plyr::mapvalues(singlerDF$Cell_subtype,
                                        from = c(as.character(metadata$Cell_subtype),"unknown"),
                                        to = c(metadata$Cell_subtype.colors,"#BEBEBE"))
table(rownames(meta.data) == rownames(singlerDF))
meta.data %<>% cbind(singlerDF)
saveRDS(meta.data, "output/Lung_time15_metadata_20220523_v2.rds")
