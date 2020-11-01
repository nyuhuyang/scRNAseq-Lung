########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#conda activate r4.0
invisible(lapply(c("Seurat","dplyr","monocle",
                   "magrittr","future","ggplot2","tidyr"), function(x) {
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
read.path = "Yang/Lung_30/Monocle2/20-ReadDE-distal.terminal.proximal.COPD-root=BC/"
cds = readRDS(paste0(read.path,basename(read.path),"_cds.rds"))
meta.data = pData(cds)
remove(cds);GC()

object = readRDS(file = "data/Lung_30_20200710.rds")
object %<>% subset(cells = rownames(meta.data))
table(rownames(object@meta.data) == rownames(meta.data))
object@meta.data = meta.data
DefaultAssay(object) = "SCT"

step <- switch(as.character(which(sapply(list(1:7, 8:12,13:29,30:77), function(x) args %in% x))),
              "1" = "1a. Determine DEGs for cells classified per state and pseudotime",
              "2" = "1b. Determine DEGs for intra-state comparisons",
              "3" = "1c. DEGs per each state (from 1 to 17) compared to other states",
              "4" = "2.Determine DEGs between groups")
print(step)
# 1a. Determine DEGs for cells classified per state and pseudotime (P+D+T+COPD combined):
if(step == "1a. Determine DEGs for cells classified per state and pseudotime"){ #  need 32 GB
    save.path = paste0(path,step,"/")
    if(!dir.exists(save.path))dir.create(save.path, recursive = T)
    
    DE_criterias = list("meta.data$State == 15 & meta.data$Pseudotime < 2",
                        "meta.data$State == 15 & meta.data$Pseudotime < 1",
                        "meta.data$State == 6",
                        "meta.data$State == 6 & meta.data$Pseudotime >13 & meta.data$Pseudotime < 14",
                        "meta.data$State == 16 & meta.data$Pseudotime >13 & meta.data$Pseudotime < 14",
                        "meta.data$State == 1 & meta.data$Pseudotime >13 & meta.data$Pseudotime < 17",
                        "meta.data$State %in% c(1,6,16,17) & meta.data$Pseudotime >13 & meta.data$Pseudotime < 14"
    )
    DE_criteria = DE_criterias[[args]]
    # write DE_criteria
    print(cat("Subset Criteria is:","\n",DE_criteria))
    cells <- rownames(meta.data)[eval(parse(text = DE_criteria))]
    object@meta.data[,"DE_criteria"] = "others"
    object@meta.data[cells,"DE_criteria"] = DE_criteria
    
    Idents(object) = "DE_criteria"
    print(table(Idents(object)))

    DEG <- FindMarkers.UMI(object, logfc.threshold = 0.01,test.use = "MAST",
                           ident.1 = DE_criteria, ident.2 = "others", 
                           only.pos = TRUE,latent.vars = "nFeature_SCT",
                           p.adjust.methods = "fdr")
    DE_criteria_name = gsub("meta.data\\$","",DE_criteria)
    DEG$cluster = DE_criteria_name
    DEG$gene = rownames(DEG)
    DEG = DEG[DEG$p_val_adj < 0.05, ]
    if(args < 10) args = paste0("0", args)
    write.csv(DEG, file = paste0(save.path,"Lung_30_monocle2.1a_",args,"_",
                                 DE_criteria_name,
                                 ".csv"))
}

# 1b. Determine DEGs for intra-state comparisons:
if(step == "1b. Determine DEGs for intra-state comparisons"){ #  need 16 GB
    args = args -7
    save.path = paste0(path,step,"/")
    if(!dir.exists(save.path))dir.create(save.path, recursive = T)
    
    DE_criterias = list(list("meta.data$State == 15 & meta.data$Pseudotime < 1",
                             "meta.data$State == 15 & meta.data$Pseudotime > 1"),
                        list("meta.data$State == 15 & meta.data$Pseudotime < 2",
                             "meta.data$State == 15 & meta.data$Pseudotime > 2"),
                        list("meta.data$State == 1 & meta.data$Pseudotime < 17",
                             "meta.data$State == 1 & meta.data$Pseudotime > 17"),
                        list("meta.data$State == 6 & meta.data$Pseudotime > 13",
                             "meta.data$State == 6 & meta.data$Pseudotime < 13"),
                        list("meta.data$State == 16 & meta.data$Pseudotime < 14",
                             "meta.data$State == 16 & meta.data$Pseudotime > 14")
    )
    DE_criteria = DE_criterias[[args]]
    # write DE_criteria
    cat("Subset Criteria is:","\n",DE_criteria[[1]],"\n",DE_criteria[[2]])
    cells_1 <- rownames(meta.data)[eval(parse(text = DE_criteria[[1]]))]
    cells_2 <- rownames(meta.data)[eval(parse(text = DE_criteria[[2]]))]
    
    object@meta.data[,"DE_criteria"] = "others"
    object@meta.data[cells_1,"DE_criteria"] = DE_criteria[[1]]
    object@meta.data[cells_2,"DE_criteria"] = DE_criteria[[2]]
    
    Idents(object) = "DE_criteria"
    object %<>% subset(idents = "others",invert = TRUE)
    print(table(Idents(object)))
    
    DEG <- FindAllMarkers.UMI(object, logfc.threshold = 0.01,test.use = "MAST",
                              only.pos = TRUE,latent.vars = "nFeature_SCT",
                              p.adjust.methods = "fdr")
    DE_criteria_name = gsub("meta.data\\$","",DE_criteria[[1]])
    
    DEG$cluster %<>% gsub("meta.data\\$","",.)
    DEG = DEG[DEG$p_val_adj < 0.05, ]
    if(args < 10) args = paste0("0", args)
    write.csv(DEG, file = paste0(save.path,"Lung_30_monocle2.1b_",args,"_",
                                 DE_criteria_name,
                                 ".csv"))
}

# 1c. DEGs per each state (from 1 to 17) compared to other states:
if(step == "1c. DEGs per each state (from 1 to 17) compared to other states"){ #  need 32 GB
    args = args -12
    save.path = paste0(path,step,"/")
    if(!dir.exists(save.path))dir.create(save.path, recursive = T)
    
    Idents(object) = "State"
    print(table(Idents(object)))
    
    DEG <- FindMarkers.UMI(object, logfc.threshold = 0.01,test.use = "MAST",
                           ident.1 = args, ident.2 = NULL, 
                           only.pos = TRUE,latent.vars = "nFeature_SCT",
                           p.adjust.methods = "fdr")
    DEG$cluster = args
    DEG$gene = rownames(DEG)
    DEG = DEG[DEG$p_val_adj < 0.05, ]
    if(args < 10) args = paste0("0", args)
    write.csv(DEG, file = paste0(save.path,"Lung_30_monocle2.1c_State=",args,".csv"))
}

# 2.Determine DEGs between groups (P vs D, P vs D+T, T vs P+D, COPD vs D) per each state (from 1 to 17).
step = "2.Determine DEGs between groups"
if(step == "2.Determine DEGs between groups"){ # need 32 GB
    save.path = paste0(path,step,"/")
    if(!dir.exists(save.path))dir.create(save.path, recursive = T)
    args = args -29
    Idents_list = list(ident1 = list("proximal",
                                     "proximal",
                                     "terminal",
                                     "COPD"),
                       ident2 = list("distal",
                                     c("distal", "terminal"),
                                     c("proximal","distal"),
                                     "distal"))
    States = c(1,5:9,12:17)
    i = ceiling((args/12) %% 4)
    if(args == 48) i = 4
    print(ident1 <- Idents_list$ident1[[i]])
    print(ident2 <- Idents_list$ident2[[i]])
    
    k = States[((args-1) %% 12)+1]
    # subset
    Idents(object) = "State"
    object %<>% subset(idents = k)
    Idents(object) = "conditions"
    object %<>% subset(idents = c(ident1,ident2))
    
    if(length(ident2) >1){
        object@meta.data$conditions = gsub(paste(ident2, collapse = "|"),
                                               paste(ident2, collapse = "+"),
                                               object@meta.data$conditions)
        ident2 = paste(ident2, collapse = "+")
    }
    Idents(object) = "conditions"
    Idents(object) %<>% factor(levels = c(ident1,ident2))
    print(table(Idents(object)))
    
    DEG <- FindAllMarkers.UMI(object, logfc.threshold = 0,test.use = "MAST",
                              only.pos = TRUE,latent.vars = "nFeature_RNA",
                              return.thresh = 1, p.adjust.methods = "fdr")
    if(args < 10) args = paste0("0", args)
    DEG = DEG[DEG$p_val_adj < 0.05, ]
    
    write.csv(DEG, file = paste0(save.path,"Lung_30_monocle2.2_",args,"_State=",k,
                                 "_",ident1,"_vs_",ident2,".csv"))
}
