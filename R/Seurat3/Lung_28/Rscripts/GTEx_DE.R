# conda activate r4.0
library(Seurat)
library(openxlsx)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
set.seed(101)

step <- switch(as.character(args),
               "1" = "1.All (both males and females) young vs old",
               "2" = "2.Males only: young vs old",
               "3" = "3.Females only: young vs old",
               "4" = "4.Young: females vs males",
               "5" = "5.Old: females vs males")
print(step)

# load data
(load(file = "data/Lung_GTEx_20200307.Rda"))
object$age =  plyr::mapvalues(object$Age.Bracket, 
                              from = c("20-29",
                                       "30-39",
                                       "40-49",
                                       "50-59",
                                       "60-69",
                                       "70-79"),
                              to = c("young",
                                     "young",
                                     "middle",
                                     "middle",
                                     "old",
                                     "old"))
Idents(object) = "age"
object %<>% subset(idents = c("young","old"))

# read optimized DEGs
DEGs <- read.csv(file = "Yang/Lung_30/DE_analysis/Optimized_cell_type_DEG.csv",row.names = 1)
cell.types <- unique(DEGs$cluster)
Idents(object) %<>% factor(levels = c("young","old"))
deg_list <- vector(mode = "list", length = length(cell.types))
names(deg_list) = cell.types

# 1.       All (both males and females):  young vs old
if(step == "1.All (both males and females) young vs old"){
    degs <- FindAllMarkers.UMI(object, 
                               #features = features,
                               logfc.threshold = 0,
                               only.pos = T,
                               return.thresh = 1, 
                               p.adjust.methods = "BH",
                               test.use = "MAST")

    for(i in 1:length(cell.types)){
        features <- DEGs[DEGs$cluster %in% cell.types[i], "gene"]
        degs1 =  degs[degs$gene %in% features & degs$cluster %in% "young", ]
        degs2 =  degs[degs$gene %in% features & degs$cluster %in% "old", ]
        degs1$p_val_adj = p.adjust(p = degs1$p_val, method = "BH",n = nrow(x = object))
        degs2$p_val_adj = p.adjust(p = degs2$p_val, method = "BH",n = nrow(x = object))
        deg_list[[i]] = rbind(degs1, degs2)
    }
}

# 2.       Males only: young vs old
if(step == "2.Males only: young vs old"){
    Idents(object) = "Sex"
    sub_object <- subset(object, idents = "male")
    Idents(sub_object) = "age"
    Idents(sub_object) %<>% factor(levels = c("young","old"))
    degs <- FindAllMarkers.UMI(sub_object, 
                               #features = features,
                               logfc.threshold = 0,
                               only.pos = T,
                               return.thresh = 1, 
                               p.adjust.methods = "BH",
                               test.use = "MAST")
    
    for(i in 1:length(cell.types)){
        features <- DEGs[DEGs$cluster %in% cell.types[i], "gene"]
        degs1 =  degs[degs$gene %in% features & degs$cluster %in% "young", ]
        degs2 =  degs[degs$gene %in% features & degs$cluster %in% "old", ]
        degs1$p_val_adj = p.adjust(p = degs1$p_val, method = "BH",n = nrow(x = sub_object))
        degs2$p_val_adj = p.adjust(p = degs2$p_val, method = "BH",n = nrow(x = sub_object))
        deg_list[[i]] = rbind(degs1, degs2)
    }
}

# 3.       Females only: young vs old
if(step == "3.Females only: young vs old"){
    Idents(object) = "Sex"
    sub_object <- subset(object, idents = "female")
    Idents(sub_object) = "age"
    Idents(sub_object) %<>% factor(levels = c("young","old"))
    degs <- FindAllMarkers.UMI(sub_object, 
                               #features = features,
                               logfc.threshold = 0,
                               only.pos = T,
                               return.thresh = 1, 
                               p.adjust.methods = "BH",
                               test.use = "MAST")
    for(i in 1:length(cell.types)){
        features <- DEGs[DEGs$cluster %in% cell.types[i], "gene"]
        degs1 =  degs[degs$gene %in% features & degs$cluster %in% "young", ]
        degs2 =  degs[degs$gene %in% features & degs$cluster %in% "old", ]
        degs1$p_val_adj = p.adjust(p = degs1$p_val, method = "BH",n = nrow(x = sub_object))
        degs2$p_val_adj = p.adjust(p = degs2$p_val, method = "BH",n = nrow(x = sub_object))
        deg_list[[i]] = rbind(degs1, degs2)
    }
}

# 4.       Young: females vs males
if(step == "4.Young: females vs males"){
    Idents(object) = "age"
    sub_object <- subset(object, idents = "young")
    Idents(sub_object) = "Sex"
    Idents(sub_object) %<>% factor(levels = c("female","male"))
    degs <- FindAllMarkers.UMI(sub_object, 
                               #features = features,
                               logfc.threshold = 0,
                               only.pos = T,
                               return.thresh = 1, 
                               p.adjust.methods = "BH",
                               test.use = "MAST")
    for(i in 1:length(cell.types)){
        features <- DEGs[DEGs$cluster %in% cell.types[i], "gene"]
        degs1 =  degs[degs$gene %in% features & degs$cluster %in% "female", ]
        degs2 =  degs[degs$gene %in% features & degs$cluster %in% "male", ]
        degs1$p_val_adj = p.adjust(p = degs1$p_val, method = "BH",n = nrow(x = sub_object))
        degs2$p_val_adj = p.adjust(p = degs2$p_val, method = "BH",n = nrow(x = sub_object))
        deg_list[[i]] = rbind(degs1, degs2)
    }
}


# 5.       Old: females vs males
if(step == "5.Old: females vs males"){
    Idents(object) = "age"
    sub_object <- subset(object, idents = "old")
    Idents(sub_object) = "Sex"
    Idents(sub_object) %<>% factor(levels = c("female","male"))
    degs <- FindAllMarkers.UMI(sub_object, 
                               #features = features,
                               logfc.threshold = 0,
                               only.pos = T,
                               return.thresh = 1, 
                               p.adjust.methods = "BH",
                               test.use = "MAST")
    for(i in 1:length(cell.types)){
        features <- DEGs[DEGs$cluster %in% cell.types[i], "gene"]
        degs1 =  degs[degs$gene %in% features & degs$cluster %in% "female", ]
        degs2 =  degs[degs$gene %in% features & degs$cluster %in% "male", ]
        degs1$p_val_adj = p.adjust(p = degs1$p_val, method = "BH",n = nrow(x = sub_object))
        degs2$p_val_adj = p.adjust(p = degs2$p_val, method = "BH",n = nrow(x = sub_object))
        deg_list[[i]] = rbind(degs1, degs2)
    }
}

write.xlsx(deg_list[lapply(deg_list,length)>0], file = paste0(path,step,".xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
