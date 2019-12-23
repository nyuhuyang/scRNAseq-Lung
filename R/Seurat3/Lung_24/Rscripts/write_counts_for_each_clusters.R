########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","magrittr","MAST"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
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

# load data
(load(file = paste0("data/Lung_24_20191206.Rda")))
object$cell.types %<>%  gsub("/","_",.)
DefaultAssay(object)  = "SCT"
Idents(object) = "cell.types"
object %<>% sortIdent
cell.types <- unique(Idents(object))
cell.type = cell.types[args]
print(paste("Expression data for =",cell.type))
sub_object <- subset(object, idents = cell.type)

# Expression data =============
meta.data = sub_object@meta.data[,c("cell.types","conditions","barcode","orig.ident")]
colnames(meta.data) %<>% sub("orig.ident","samples",.)
colnames(meta.data) %<>% sub("conditions","regions",.)
rownames(meta.data) = meta.data$barcode
data = DelayedArray::t(sub_object@assays$SCT@data)
data = merge(meta.data,data,by="row.names",all.x=TRUE)
data %<>% arrange(match(regions, c("proximal","distal","terminal")))
rownames(data) = data$Row.names
data = data[,-grep("Row.names|cell.types|barcode",colnames(data))]
write.csv(DelayedArray::t(data), paste0(path,"Lung_24-",args,"_",cell.type,".csv"))


# test expression data
#exp <- read.csv("Yang/proximal_distal_terminal/Non-Integration/Counts and meta/Cell types/Lung_24-12_Basal cells.csv",row.names = 1, stringsAsFactors = F)
#genes <- c("S100A9","ALDH3A1","FST","PLAU","KRT14","MMP10","SFN",
#           "S100A8","GJB2","LYZ","TNC","KRT6A","AKR1C2","DKK1","S100A2",
#           "PRR4","S100A11","FGFBP1")
#regions <- exp["regions",]  %>% unlist
#table(regions)
#ident.1 = names(regions)[regions %in% "proximal"]
#ident.2 = names(regions)[regions %in% c("distal","terminal")]
#exp_genes <- data.matrix(exp[genes,])
#de.results <- FindMarkers(
#        object = exp_genes,
#        cells.1 = ident.1,
#        cells.2 = ident.2,
#        features = genes,
#        test.use = "MAST")
