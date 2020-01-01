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

# select conditions
conditions = c("proximal","distal","terminal")
(con <- conditions[args])
# read file
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
# Load Seurat
(load(file = paste0("data/Epi_24_",con,"_20191223.Rda")))
object$cell_types <- plyr::mapvalues(object$cell_types,
                                     from = df_cell_types$`Cell types`,
                                     to = df_cell_types$Abbreviation)
DefaultAssay(object)  = "SCT"

# Expression data =============
meta.data = object@meta.data[,c("cell_types","conditions","barcode","orig.ident")]
colnames(meta.data) %<>% sub("orig.ident","samples",.)
colnames(meta.data) %<>% sub("conditions","regions",.)
rownames(meta.data) = meta.data$barcode
data = DelayedArray::t(object@assays$SCT@data)
data = merge(meta.data,data,by="row.names",all.x=TRUE)
data %<>% arrange(match(regions, c("proximal","distal","terminal")))
rownames(data) = data$Row.names
data = data[,-grep("Row.names|barcode",colnames(data))]
write.csv(DelayedArray::t(data), paste0(path,"Lung_24-",args,"_",cell.type,".csv"))