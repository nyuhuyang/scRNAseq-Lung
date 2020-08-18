########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load files
object = readRDS(file = "data/Lung_30_20200710.rds") 
# Need 64GB
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")

step = "A - Sample types"

if(step == "A - Sample types"){ # need 32 GB
    Idents_list = list(ident1 = list("proximal",
                                     "proximal",
                                     "distal",
                                     "distal",
                                     "terminal",
                                     "COPD",
                                     "COPD"),
                       ident2 = list("distal",
                                     c("distal", "terminal"),
                                     c("proximal","terminal"),
                                     "terminal",
                                     c("proximal","distal"),
                                     "distal",
                                     c("distal","terminal")))
    i = ceiling((args/54) %% 7)
    if(args == 378) i = 7
    print(ident1 <- Idents_list$ident1[[i]])
    print(ident2 <- Idents_list$ident2[[i]])
    
    k = ((args-1) %% 54)+1
    
    df_annotation <- readxl::read_excel("doc/20200815_Comparison groups of cells - plan.xlsx",
                                        sheet = step)
    groups = df_annotation$`CELL GROUPS`
    groups = groups[!is.na(groups)]
    groups = gsub(" \\(.*", "", groups)
    group_list <- stringr::str_split(groups, pattern = "\\+")
    names(group_list) = groups
    group_list[group_list == "ALL IMMUNE CELLS"][[1]] = 
        unique(unlist(group_list[35:53]))
    group_list[group_list == "ALL CELLS"][[1]] = 
        unique(unlist(group_list[2:53]))
    print(group <- group_list[[k]])
    
    Idents(object)= "annotations3"
    sub_object <- subset(object, idents = group)
    rm(object);GC()
    Idents(sub_object) = "conditions"
    sub_object <- subset(sub_object, idents = c(ident1,ident2))
    if(length(ident2) >1){
        sub_object@meta.data$conditions = gsub(paste(ident2, collapse = "|"),
                                               paste(ident2, collapse = "+"),
                                               sub_object@meta.data$conditions)

    }
    ident2 = paste(ident2, collapse = "+")
    Idents(sub_object) = "conditions"
    Idents(sub_object) %<>% factor(levels = c(ident1,ident2))
    DEG <- FindAllMarkers.UMI(sub_object, logfc.threshold = 0,test.use = "MAST",
                              return.thresh = 0.05, p.adjust.methods = "fdr")
    write.csv(DEG, file = paste0(path,"Lung_30_",args,"_celltypes=",k,
                                 "_",ident1,"_vs_",ident2,".csv"))
}
