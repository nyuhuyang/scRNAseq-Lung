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
DefaultAssay(object)  = "SCT"
Idents(object) = "cell.types"
object %<>% sortIdent
cell.types <- unique(Idents(object))
cell.type = cell.types[args]
print(paste("Expression data for =",cell.type))
sub_object <- subset(object, idents = cell.type)

# Expression data =============
meta.data = sub_object@meta.data[,c("cell.types","conditions","barcode")]
colnames(meta.data) %<>% sub("conditions","regions",.)
meta.data %<>% arrange(match(regions, c("proximal","distal","terminal")))
rownames(meta.data) = meta.data$barcode
data = DelayedArray::t(sub_object@assays$SCT@data)
data = cbind(meta.data,data)
data = data[,-grep("cell.types|barcode",colnames(data))]
write.table(DelayedArray::t(data), paste0(path,"Lung_24-",args,"_",
                                          cell.type,"_counts.txt"),
            sep='\t', quote=F)