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
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

cell.types = c("AT1","AT2-a","AT2-b","B","B1","BC",     
                "BC-p","C1",  "C2",  "C3",  "C4",  "Cr", 
                "DC-p","DC1", "DC2", "En-A","En-C","En-C1",
                "En-L","En-SM","En-V","F1",  "F2",  "F3", 
                "F4",  "F5",  "Gli", "H",   "IC",  "Ion",
                "M-p", "M0",  "M1",  "M2",  "MC",  "MEC",
                "MEC-1","Mon-1","Mon-2","NEC", "Neu-1","Neu-2",
                "Nr",  "p-C", "PC",  "Pr",  "S",   "S-d",
                "SM1", "SM2", "SM3", "SMG-Muc","SMG-Ser","T-cn",
                "T-inf","T-NK","T-p", "T-reg","T-rm","Un-AT-cont",
                "Un-En-cont","Un-F-cont","Un-M-cont","Un-Prox-6", "Un-S-cont","Un-SAE-cont",
                "Un-SM-cont","Un-T-cont") 
(c = cell.types[args])
doc.path <- paste0("Yang/proximal_distal_terminal_COPD/Harmony/Annotations/")
if(!dir.exists(doc.path))dir.create(doc.path, recursive = T)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#======== read csv =================
(csv_files <-list.files(doc.path, pattern = "20200206.csv"))
anno_files <- lapply(paste0(doc.path, csv_files), read.csv, stringsAsFactors = F)
anno <- rbindlist(anno_files) %>% as.data.frame()
#======== adjust cell labels =================
anno = anno[order(anno$cell.labels),];dim(anno)
df1 = anno[!duplicated(anno$barcodes),];dim(df1)
df2 = anno[duplicated(anno$barcodes),];dim(df2)

df2$samples = gsub("_.*", "", df2$barcodes)
df2 = df2[,c("barcodes", "samples", "cell.labels")]
df_ambigous  = inner_join(df1, df2, by ="barcodes")
# ambigous have two labels
df_ambigous %<>% .[.[,"cell.labels.x"] != .[,"cell.labels.y"],];dim(df_ambigous)

df_unique <- anno[!anno$barcodes %in% df_ambigous$barcodes,];dim(df_unique)
# unique have one label or two same labels
df_unique <- df_unique[!duplicated(df_unique$barcodes),];dim(df_unique)
# unique contains known and unkown
unknown <- grepl("^Un-",df_unique$cell.labels);table(unknown)
df_unknown = df_unique[unknown,];dim(df_unknown)
df_known = df_unique[!unknown,];dim(df_known)

object = readRDS(file = "data/Lung_28_Global_20200206.rds") 
#genes = readxl::read_excel("doc/Genes-for-ambiguous-cell-annotation.xlsx", 
#                           col_names = F) %>% pull
#exp = object[["SCT"]]@data[genes,df_ambigous$barcodes]
#df_ambigous %<>% cbind(Matrix::t(exp))
meta.data = object@meta.data;dim(meta.data)
meta.data$barcodes = rownames(meta.data)

df_nolabel = meta.data[!meta.data$barcodes %in% anno$barcodes,];dim(df_nolabel)
ambigous_unkown_nolabel = unique(c(df_ambigous$barcodes,
                                   df_unknown$barcodes,
                                   df_nolabel$barcodes))
length(ambigous_unkown_nolabel)
#write.csv(ambigous_unkown_nolabel, paste0(doc.path,"Lung_28_ambigous_unkown_nolabel_barcodes.csv"))

meta.data = inner_join(meta.data, df_unique, by = "barcodes");dim(meta.data)
rownames(meta.data) = meta.data$barcodes
object %<>% subset(cells = rownames(meta.data))
table(rownames(meta.data) == colnames(object))
object[["cell.labels"]] = meta.data$cell.labels
Idents(object) = "cell.labels"
object %<>% sortIdent()
object %<>% AddMetaColor(label= "cell.labels",colors = Singler.colors)
UMAPPlot.1(object, label = T, cols = Singler.colors,
           label.repel = T, no.legend = T,
           do.print = T, do.return = F)

system.time(Lung_markers <- FindMarkers.UMI(object, ident.1 = c, ident.2 = NULL,
                                               logfc.threshold = 0.5, only.pos = T,
                                               test.use = "MAST",min.cells.group = 2))
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
if(args < 10) args = paste0("0", args)
write.csv(Lung_markers,paste0(path,"Lung_28-",args,"_",c,".csv"))
