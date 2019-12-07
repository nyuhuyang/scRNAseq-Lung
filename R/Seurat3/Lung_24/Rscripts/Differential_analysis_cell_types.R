########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","tidyr","magrittr","gplots","MAST"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

df_labels_Pairs <- readxl::read_excel("doc/DEG cell types comparison plan.xlsx")
df_labels_Pairs = df_labels_Pairs[!is.na(df_labels_Pairs$`Group 1`),]
labels_Pairs = df_labels_Pairs[args,]
(ident1 <- strsplit(labels_Pairs$`Group 1`, "\\+") %>% unlist())
(ident2 <- strsplit(labels_Pairs$`Group 2`, "\\+")%>% unlist())


# load data
(load(file = "data/Lung_24_20191206.Rda"))
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
table(df_cell_types$`Cell types` %in% object$cell.types)

object$cell.types %<>% plyr::mapvalues(from = df_cell_types$`Cell types`,
                                       to = df_cell_types$Abbreviation)
DefaultAssay(object) = "RNA"
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"
Lung_markers <- FindMarkers.UMI(object,
                                ident.1 = ident1, 
                                ident.2 = ident2,
                                logfc.threshold = 0, only.pos = F,
                                test.use = "MAST")

write.csv(Lung_markers,paste0(path,"Lung_24_RNA_pairwise_",args,"_",
                              paste(ident1,collapse = "."),
                              ".csv"))

