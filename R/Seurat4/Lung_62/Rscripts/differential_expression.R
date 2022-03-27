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

object = readRDS(file ="data/Lung_62_20220322.rds")

(step = c("resolutions","Cell_subtype")[2])

if(step == "resolutions"){# 32GB
    opts = data.frame(ident = c(rep("SCT_snn_res.0.01",33),
                                rep("SCT_snn_res.0.1",55),
                                rep("SCT_snn_res.0.2",67),
                                rep("SCT_snn_res.0.5",99),
                                rep("SCT_snn_res.0.8",122)),
                      num = c(0:32,
                              0:54,
                              0:66,
                              0:98,
                              0:121)
                      )

    opt = opts[args,]
    print(opt)
    #==========================
    Idents(object) = opt$ident

    markers = FindMarkers_UMI(object, ident.1 = opt$num,
                              group.by = opt$ident,
                              assay = "SCT",
                              #min.pct = 0.01,
                              logfc.threshold = 0.5,
                                 only.pos = T#,
                                 #test.use = "MAST",
                                 #latent.vars = "nFeature_SCT"
                              )
    markers$cluster = opt$num
    num = opt$num
    if(args < 10) num = paste0("0",num)
    if(args < 100) num = paste0("0",num)

    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)

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
