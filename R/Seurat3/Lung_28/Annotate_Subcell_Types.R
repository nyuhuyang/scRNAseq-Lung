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
df_ambigous = inner_join(df2, df1, by = "barcodes");dim(df_ambigous)
df_ambigous %<>% .[.[,"cell.labels.x"] != .[,"cell.labels.y"],];dim(df_ambigous)

object = readRDS(file = "data/Lung_28_Global_20200206.rds") 
genes = readxl::read_excel("doc/Genes-for-ambiguous-cell-annotation.xlsx", 
                           col_names = F) %>% pull
exp = object[["SCT"]]@data[genes,df_ambigous$barcodes]
df_ambigous %<>% cbind(Matrix::t(exp))

df_ambigous_old <- read.csv(paste0(path,"Lung_28_amb_annotations_old.csv"),
                            stringsAsFactors = F,row.names = 1)
dim(df_ambigous_old)
dim(df_ambigous)

table(df_ambigous$barcodes %in% df_ambigous_old$barcodes)
df_ambigous %<>% .[!(.[,"barcodes"] %in% df_ambigous_old[,"barcodes"]),]

write.csv(df_ambigous, paste0(path,"Lung_28_amb_annotations_new.csv"))


df_full = full_join(df1, df2, by = "barcodes");dim(df_full)
meta.data = object@meta.data
meta.data$barcodes = rownames(meta.data)
meta.data = full_join(meta.data, df_full, by = "barcodes")
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

df_nolabel = meta.data[is.na(meta.data$cell.labels.x),c("barcodes","orig.ident")]
exp = object[["SCT"]]@data[genes,df_nolabel$barcodes]
df_nolabel %<>% cbind(Matrix::t(exp))
write.csv(df_nolabel, paste0(path,"Lung_28_no_annotations.csv"))

ambigous_nolabel = unique(c(df_ambigous_old$barcodes, 
                            df_ambigous$barcodes, 
                            df_nolabel$barcodes))
write.csv(ambigous_nolabel, paste0(path,"Lung_28_ambigous_nolabel_barcodes.csv"))
