########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","magrittr","MAST"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/Expression/")
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
object$cell_types <- plyr::mapvalues(object$cell.types,
                                     from = df_cell_types$`Cell types`,
                                     to = df_cell_types$Abbreviation)
DefaultAssay(object)  = "SCT"
data = object@assays$SCT@data
# 1. gene list =============
write.csv(rownames(data), paste0(path,"Epi_24-",con,"_gene_list.csv"),
          row.names = FALSE,col.names = FALSE)

# 2. Cell grouping =========
meta.data = object@meta.data[,c("cell_types","orig.ident")]
colnames(meta.data) %<>% sub("orig.ident","samples",.)
table(rownames(meta.data) == colnames(data))
rownames(meta.data) = NULL
write.csv(t(meta.data), paste0(path,"Epi_24-",con,"_Cell_grouping.csv"),
          row.names = TRUE,col.names = FALSE)
# Expression data =============
rownames(data) =NULL
colnames(data) =NULL
write.csv(data, paste0(path,"Epi_24-",con,"_exp.csv"),
          row.names = FALSE,col.names = FALSE)


#data = DelayedArray::t(object@assays$SCT@data)
#data = merge(meta.data,data,by="row.names",all.x=TRUE)
#data %<>% arrange(match(regions, c("proximal","distal","terminal")))
#rownames(data) = data$Row.names
#data = data[,-grep("Row.names|barcode",colnames(data))]
#write.csv(DelayedArray::t(data), paste0(path,"Epi_24-",con,"_exp.csv"))