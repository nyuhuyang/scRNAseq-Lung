library(Seurat)
library(magrittr)
library(dplyr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

object = readRDS(file ="data/Lung_SCT_63_20220408.rds")
object@meta.data = readRDS(file = "data/Lung_63_20220408_meta.data_v2.rds")
(step = c("resolutions","Cell_subtype")[1])

if(step == "resolutions"){# 32GB
    opts = data.frame(ident = c(rep("SCT_snn_res.0.01",6),
                                rep("SCT_snn_res.0.1",22),
                                rep("SCT_snn_res.0.2",30),
                                rep("SCT_snn_res.0.5",48),
                                rep("SCT_snn_res.0.8",66),#
                                rep("SCT_snn_res.0.9",68),
                                rep("SCT_snn_res.1",72),
                                rep("SCT_snn_res.2",107),
                                rep("SCT_snn_res.3",129),
                                rep("SCT_snn_res.4",151),
                                rep("SCT_snn_res.5",180)),
                      num = c(0:5,
                              0:21,
                              0:29,
                              0:47,
                              0:65,
                              0:67,
                              0:71,
                              0:106,
                              0:128,
                              0:150,
                              0:179)
                      )

    opt = opts[args,]
    print(opt)
    #==========================
    Idents(object) = opt$ident

    markers = FindMarkers_UMI(object, ident.1 = opt$num,
                              group.by = opt$ident,
                              assay = "SCT",
                              #min.pct = 0.01,
                              logfc.threshold = 0.25,
                                 only.pos = T#,
                                 #test.use = "MAST",
                                 #latent.vars = "nFeature_SCT"
                              )
    markers$cluster = opt$num
    num = opt$num
    if(args < 10) num = paste0("0",num)
    if(args < 100) num = paste0("0",num)
    if(args < 1000) num = paste0("0",num)
    
    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)
    if(args < 1000) arg = paste0("0",arg)
    
    write.csv(markers,paste0(path,arg,"_",opt$ident,"_",num, ".csv"))
}


if(step == "Cell_subtype"){# 32~64GB
   DefaultAssay(object) = "SCT"
   object %<>% subset(subset = Cell_subtype != "Un" &
                           orig.ident %in% c("WC_30_T_Adj_Dex",
                                             "WC_30_T_Adj_Cont")
    )
   Cell_subtypes = c("All",sort(unique(object$Cell_subtype)))
   type = Cell_subtypes[args]
    if(type != "All") object %<>% subset(subset = Cell_subtype == type)
    #==========================
    Idents(object) = "orig.ident"
    print(type)
    
    markers = FindMarkers_UMI(object, ident.1 = "WC_30_T_Adj_Dex",
                              ident.2 = "WC_30_T_Adj_Cont",
                              group.by = "orig.ident",
                              assay = "SCT",
                              min.pct = 0.01,
                              logfc.threshold = 0.1,
                                 only.pos = F)
    markers$gene = rownames(markers)
    markers$cluster = type
    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)

    save_path <- paste0(path,step,"/")
    if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
    
    write.csv(markers,paste0(save_path,arg,"-",type,"-Dex_vs_Count",".csv"))
}
