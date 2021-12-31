# conda activate r4.0.3
library(Seurat)
library(magrittr)
library(kableExtra)
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
# create singleR data frame
###############################
pred = readRDS("output/Lung_SCT_20211201_singleR_bulk_krasnowLung.rds")
pred = readRDS("output/Lung_SCT_20211201_singleR_sc_krasnowLung.rds")

meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")


singlerDF = data.frame("azimuth_annotation" = pred$pruned.labels,
                       row.names = rownames(pred))
table(is.na(singlerDF$azimuth_annotation))
singlerDF$azimuth_annotation[is.na(singlerDF$azimuth_annotation)]= "unknown"
table(singlerDF$azimuth_annotation)
##############################
# adjust cell label
##############################
# combine cell types
path = "../seurat_resources/azimuth/human_lung/"
metadata = fread(paste0(path,"krasnow_hlca_10x_metadata.csv.gz"))
metadata = metadata[!duplicated(metadata$free_annotation),]

singlerDF$compartment = plyr::mapvalues(singlerDF$azimuth_annotation,
                                        from = metadata$free_annotation,
                                        to = metadata$compartment)
table(rownames(meta.data) == rownames(singlerDF))
meta.data %<>% cbind(singlerDF)
saveRDS(meta.data, "output/20211209/meta.data_azimuth_subtype.rds")
