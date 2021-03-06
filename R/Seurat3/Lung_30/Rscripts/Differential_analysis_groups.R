########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","stringi",
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

# load files
object = readRDS(file = "data/Lung_30_20200710.rds") 
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")

step = "surface_airway_epithelial"
save.path <- paste0("Yang/Lung_30/DE_analysis/",step,"/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

if(step == "A_Sample_types"){ #  need 64 GB for all cells. need 32 GB for others.
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


if(step == "B_Cell_groups"){ # need 32 GB
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


if(step == "A_Sample_types-age"){ #  need 64 GB for all cells. need 32 GB for others.
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

if(step == "F_EVGs"){ # need 64 GB for all cells. need 32 GB for others.
    cell.types = c("AT1","AT2","AT2-1","AT2-p","BC-p","BC1","BC2","C1","C2","C3",
                   "d-S","g-Muc","g-Ser","H","IC1","IC2","IC3","Ion","ME","NE",
                   "p-C","S","Cr","En-a","En-c","En-c1","En-l","En-p","En-sm","En-v",
                   "F1","F2","F3","F4","Gli","Nr","Pr","SM1","SM2","SM3",
                   "B","DC","M-p","M0","M1","M1-2","M2","MC","Mon","Neu",
                   "p-DC","PC","T-cn","T-ifn","T-NK","T-p","T-reg","T-rm","T7")
    len_1 = length(cell.types)
    Idents_list = list(ident1 = list("proximal",
                                     "distal",
                                     "distal",
                                     "terminal",
                                     "COPD",
                                     "older" = c("UNC-54-D", "UNC-57-D", "UNC-66-D", "UNC-70-D")),
                       ident2 = list(c("distal", "terminal"),
                                     "proximal",
                                     c("proximal","terminal"),
                                     c("proximal","distal"),
                                     "distal",
                                     "younger" = c("UNC-44-D", "UNC-48-D", "UNC-55-D", "UNC-67-D", "UNC-69-D", "UNC-71-D", "VU-27-D")))
    len_2 = length(Idents_list$ident1)
    
    i = ceiling((args/len_1) %% len_2)
    if(args == len_1*len_2) i = len_2 #max(len_1*len_2)=354
    
    print(ident1 <- Idents_list$ident1[[i]])
    print(ident2 <- Idents_list$ident2[[i]])
    
    print(cell.type <-cell.types[((args-1) %% len_1)+1])

    # change cell name
    anno <- readxl::read_excel("doc/Annotations/Cell type abbreviation.xlsx")
    Idents(object) = "annotations3"
    object$cell.types = plyr::mapvalues(object$annotations3,
                                        from = anno$Abbreviation,
                                        to = anno$`Revised abbreviations`)
    
    # subset
    compare = ifelse(i<6,"conditions","age")
    
    if(i==6) {
        object$age = ""
        object@meta.data[object$orig.ident %in% ident1,"age"] = "older"
        object@meta.data[object$orig.ident %in% ident2,"age"] = "younger"
        ident1="older"
        ident2="younger"
    }
    Idents(object) = compare
    object <- subset(object, idents = c(ident1,ident2))
    
    Idents(object) = "cell.types"
    object %<>% subset(idents = cell.type)


    if(length(ident2) >1){
        object@meta.data[,compare] = gsub(paste(ident2, collapse = "|"),
                                          paste(ident2, collapse = "+"),
                                          object@meta.data[,compare])
        ident2 = paste(ident2, collapse = "+")
    }
    Idents(object) = compare
    Idents(object) %<>% factor(levels = c(ident1,ident2))
    DEG <- FindAllMarkers.UMI(object, logfc.threshold = 0,test.use = "MAST",
                              only.pos = FALSE,latent.vars = "nFeature_RNA",
                              return.thresh = 1, p.adjust.methods = "fdr")
    if(args < 10) args = paste0("0", args)
    write.csv(DEG, file = paste0(save.path,"Lung_30_",args,"_",i,
                                 "_",cell.type,"_",ident1,"_vs_",ident2,".csv"))
}

if(step == "groups"){
    anno <- readxl::read_excel("doc/Annotations/Cell type abbreviation.xlsx")
    table(anno$Abbreviation %in% object$annotations3)
    object$cell_types <- plyr::mapvalues(object$annotations3,
                                         from = anno$Abbreviation,
                                         to = anno$`Revised abbreviations`)
    cell.type_list <- list("Epithelial" = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","d-S",
                                            "H","p-C","C1","C2","C3","Ion","NE","g-Muc",
                                            "g-Ser","AT1","AT2","AT2-1","AT2-p","ME"),
                           "Surface Airway Epithelial" = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","d-S",
                                                       "H","p-C","C1","C2","C3","Ion","NE"),
                           "Airway Epithelial" = c("AT1","AT2","AT2-1","AT2-p"),
                           "Structural" = c("ME","Cr","Gli","F1","F2","F3","F4",
                                            "Nr","Pr","SM1","SM2","SM3",
                                            "En-a","En-c","En-c1","En-l","En-p","En-sm","En-v"),
                           "Immune" = c("B","DC",
                                        "M-p","M0","M1","M1-2","M2","MC","Mon","Neu",
                                        "p-DC","PC","RBC",
                                        "T-cn","T-ifn","T-int","T-NK","T-p","T-reg","T-rm")
                           )
    
    cell.type_unlist = unlist(cell.type_list)
    cell.type = cell.type_unlist[args]
    group = gsub("[0-9+]","",names(cell.type))
    
    Idents(object) = "cell_types"
    object %<>% subset(idents = cell.type_list[[group]])
    GC()
    Lung_markers <- FindMarkers.UMI(object, ident.1 = cell.type, ident.2 = NULL,
                                    logfc.threshold = 0, only.pos = F,
                                    latent.vars = "nFeature_RNA",
                                    test.use = "MAST",min.cells.group = 2)
    Lung_markers$gene = rownames(Lung_markers)
    Lung_markers$cluster = paste(cell.type, "vs.", group)
    if(args < 10) args = paste0("0", args)
    
    save.path <- paste0(save.path,group,"/")
    if(!dir.exists(save.path))dir.create(save.path, recursive = T)
    write.csv(Lung_markers,paste0(save.path,"Lung_30-",args,"_FC0_",cell.type,"_vs_", group,".csv"))
}

if(step == "surface_airway_epithelial"){ 
    object$cell_types %<>% gsub("d-S","TASC",.)
    cell.type_list <- list("Surface Airway Epithelial" = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","TASC",
                                                           "H","p-C","C1","C2","C3","Ion","NE")
    )
    Idents(object) = "cell_types"
    object %<>% subset(idents = cell.type_list$`Surface Airway Epithelial`)
    Idents(object) = "conditions"
    Idents_list = list(ident1 = list("distal",
                                     c("distal","terminal"),
                                     "COPD"),
                       ident2 = list("proximal",
                                     "proximal",
                                     "distal"))
    print(ident1 <- Idents_list$ident1[[args]])
    print(ident2 <- Idents_list$ident2[[args]])
    
    GC()
    DEG <- FindMarkers.UMI(object, ident.1 = ident1, ident.2 = ident2,features = c("SCGB3A2","SCGB1A1"),
                                    logfc.threshold = 0, only.pos = F, min.pct = 0, # need 64GB
                                    latent.vars = "nFeature_RNA",
                                    test.use = "MAST")
    if(length(ident1) >1) ident1 %<>% paste(collapse = "+")
    DEG$cluster = paste(ident1, "vs.", ident2)
    DEG$gene = rownames(DEG)
    DEG = DEG[order(DEG$avg_logFC,decreasing = T),]
    if(args < 10) args = paste0("0", args)
    write.csv(DEG, file = paste0(save.path,"Lung_30_",args,
                                 "_",ident1,"_vs_",ident2,".csv"))
    }
