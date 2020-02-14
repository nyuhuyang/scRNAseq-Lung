########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("Yang/proximal_distal_terminal_COPD/Harmony/Annotations/")
if(!dir.exists(path))dir.create(path, recursive = T)

#======== read csv =================
(csv_files <-list.files(path, pattern = "20200206.csv"))
anno_files <- lapply(paste0(path, csv_files), read.csv, stringsAsFactors = F)
anno <- rbindlist(anno_files) %>% as.data.frame()
#======== adjust cell labels =================
anno = anno[order(anno$cell.labels),]
df1 = anno[!duplicated(anno$barcodes),]
df2 = anno[duplicated(anno$barcodes),]

df2$samples = gsub("_.*", "", df2$barcodes)
df2 = df2[,c("barcodes", "samples", "cell.labels")]
df_s = left_join(df2, df1, by = "barcodes")

object = readRDS(file = "data/Lung_28_Global_20200206.rds") 
#<- readxl::read_excel("doc/20190815_scRNAseq_info.xlsx")
genes = readxl::read_excel("doc/Genes-for-ambiguous-cell-annotation.xlsx", 
                           col_names = F) %>% pull
exp = object[["SCT"]]@data[genes,df_s$barcodes]
df_s %<>% cbind(Matrix::t(exp))
write.csv(df_s, paste0(path,"Lung_28_amb_annotations.csv"))

df = full_join(df1, df2, by = "barcodes")
meta.data = object@meta.data
meta.data$barcodes = rownames(meta.data)
meta.data = full_join(meta.data, df, by = "barcodes")
meta.data = meta.data[!duplicated(meta.data$barcodes),]
rownames(meta.data) = meta.data$barcodes
meta.data = meta.data[colnames(object),]
object[["cell.labels"]] = meta.data$cell.labels.x
Idents(object) = "cell.labels"
object %<>% sortIdent()
object %<>% AddMetaColor(label= "cell.labels",colors = Singler.colors)
UMAPPlot.1(object, label = T, cols = Singler.colors,
           label.repel = T, no.legend = T,
           do.print = T, do.return = F)
