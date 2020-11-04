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
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")

step = "B - Cell groups"

if(step == "A - Sample types"){ #  need 64 GB for all cells. need 32 GB for others.
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
    i = ceiling((args/56) %% 7)
    if(args == 392) i = 7
    print(ident1 <- Idents_list$ident1[[i]])
    print(ident2 <- Idents_list$ident2[[i]])
    
    k = ((args-1) %% 56)+1
    
    df_annotation <- readxl::read_excel("doc/20200903_Comparison groups of cells - Yang modified.xlsx",
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
    object <- subset(object, idents = group)
    Idents(object) = "conditions"
    object <- subset(object, idents = c(ident1,ident2))
    if(length(ident2) >1){
        object@meta.data$conditions = gsub(paste(ident2, collapse = "|"),
                                               paste(ident2, collapse = "+"),
                                               object@meta.data$conditions)
        ident2 = paste(ident2, collapse = "+")
    }
    Idents(object) = "conditions"
    Idents(object) %<>% factor(levels = c(ident1,ident2))
    DEG <- FindAllMarkers.UMI(object, logfc.threshold = 0,test.use = "MAST",
                              only.pos = FALSE,
                              return.thresh = 1, p.adjust.methods = "fdr")
    if(args < 10) args = paste0("0", args)
    write.csv(DEG, file = paste0(path,"Lung_30_A_",args,"_celltypes=",k,
                                 "_",groups[k],"_",ident1,"_vs_",ident2,".csv"))
}


if(step == "B - Cell groups"){ # need 32 GB
    df_annotation <- readxl::read_excel("doc/20200903_Comparison groups of cells - Yang modified.xlsx",
                                        sheet = step)
    df_annotation = df_annotation[!is.na(df_annotation$`Group 1`),]

    groups = df_annotation[args,]
    groups$`Group 1` %<>% gsub(" \\(.*\\)", "", .) %>% gsub(" ", "", .)
    groups$`Group 2` %<>% gsub(" \\(.*\\)", "", .) %>% gsub(" ", "", .)
    print(ident1 <- stringr::str_split(groups$`Group 1`, pattern = "\\+")[[1]])
    print(ident2 <- stringr::str_split(groups$`Group 2`, pattern = "\\+")[[1]])
    
    Idents(object)= "annotations3"
    object <- subset(object, idents = c(ident1,ident2))

    object$annotations3 = gsub(paste(ident1, collapse = "|"),
                               paste(ident1, collapse = "+"),
                               object$annotations3)
    ident1 = paste(ident1, collapse = "+")

    object$annotations3 = gsub(paste(ident2, collapse = "|"),
                               paste(ident2, collapse = "+"),
                               object$annotations3)
    ident2 = paste(ident2, collapse = "+")
    
    Idents(object) = "annotations3"
    Idents(object) %<>% factor(levels = c(ident1,ident2))
    DEG <- FindAllMarkers.UMI(object, logfc.threshold = 0,test.use = "MAST",
                              return.thresh = 0.05, only.pos = FALSE,
                              p.adjust.methods = "fdr")
    if(args < 10) args = paste0("0", args)
    write.csv(DEG, file = paste0(path,"Lung_30_B_",args,"_",
                                 ident1,"_vs_",ident2,".csv"))
    }


if(step == "A - Sample types-age"){ #  need 64 GB for all cells. need 32 GB for others.
    Idents_list = list("older" = c("UNC-54-D", "UNC-57-D", "UNC-66-D", "UNC-70-D"),
                       "younger" = c("UNC-44-D", "UNC-48-D", "UNC-55-D", "UNC-67-D", "UNC-69-D", "UNC-71-D", "VU-27-D"))
    Idents(object) = "orig.ident"
    object %<>% subset(idents = unlist(Idents_list))
    object@meta.data$age = "younger"
    object@meta.data[object$orig.ident %in% Idents_list$older,"age"] = "older"

    df_annotation <- readxl::read_excel("doc/20200903_Comparison groups of cells - Yang modified.xlsx",
                                        sheet = sub("-age","",step))
    groups = df_annotation$`CELL GROUPS`
    groups = groups[!is.na(groups)]
    groups = gsub(" \\(.*", "", groups)
    group_list <- stringr::str_split(groups, pattern = "\\+")
    names(group_list) = groups
    group_list[group_list == "ALL IMMUNE CELLS"][[1]] = 
        unique(unlist(group_list[35:53]))
    group_list[group_list == "ALL CELLS"][[1]] = 
        unique(unlist(group_list[2:53]))
    print(group <- group_list[[args]])
    
    Idents(object)= "annotations3"
    object <- subset(object, idents = group)

    Idents(object) = "age"
    Idents(object) %<>% factor(levels = c("older","younger"))
    DEG <- FindAllMarkers.UMI(object, logfc.threshold = 0,test.use = "MAST",
                              only.pos = FALSE,
                              return.thresh = 1, p.adjust.methods = "fdr")
    if(args < 10) args = paste0("0", args)
    write.csv(DEG, file = paste0(path,"Lung_30_A_age_",args,"_celltypes=",groups[args],
                                 ".csv"))
}