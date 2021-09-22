library(Seurat)
library(magrittr)
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

object = readRDS(file = "data/Lung_SCT_30_20210831.rds")

step = c("resolutions","cell_types")[2]
if(step == "resolutions"){
    opts = data.frame(ident = c(rep("UMAP_res=0.8",84),
                                rep("UMAP_res=1",93),
                                rep("UMAP_res=2",137),
                                rep("UMAP_res=3",174),
                                rep("UMAP_res=4",200),
                                rep("UMAP_res=5",227)),
                      num = c(0:83,
                              0:92,
                              0:136,
                              0:173,
                              0:199,
                              0:226)
                      )

    opt = opts[args,]
    print(opt)
    #==========================
    meta.data = object@meta.data[,grep("SCT_snn_res",colnames(object@meta.data),invert = TRUE)]
    object@meta.data = meta.data

    spread <- 1.4
    min.dist <- 0.5
    file.name = paste0("dist.",min.dist,"_spread.",spread)

    umap = readRDS(paste0("output/20210901/","umap_",file.name,".rds"))
    meta.data = readRDS(paste0("output/20210901/meta.data_",file.name,".rds"))
    colnames(meta.data) %<>% gsub(paste0(file.name,"."),"UMAP_res=",.)

    object@meta.data %<>% cbind(meta.data)
    Idents(object) = opt$ident

    markers = FindMarkers_UMI(object, ident.1 = opt$num,
                              group.by = opt$ident,
                              assay = "SCT",
                              logfc.threshold = 0.5,
                                 only.pos = T,
                                 test.use = "MAST",
                                 latent.vars = "nFeature_SCT")
    markers$cluster = opt$num
    num = opt$num
    if(args < 10) num = paste0("0",num)
    if(args < 100) num = paste0("0",num)

    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)

    write.csv(markers,paste0(path,arg,"_",opt$ident,"_",num, ".csv"))
}

if(step == "cell_types"){
    opts = data.frame(ident = c(rep("Cell_subtype",47),
                                rep("Cell_type",31),
                                rep("UMAP_land",20),
                                rep("Family",7),
                                rep("Superfamily",3)),
                      num = c(1:47,
                              1:31,
                              1:20,
                              1:7,
                              1:3)
                      )

    opt = opts[args,]

    #==========================
    object %<>% subset(subset = Cell_subtype != "Un"
                       &  Doublets == "Singlet"
    )
    Idents(object) = opt$ident
    opt$type = sort(levels(object))[opt$num]
    print(opt)
    

    markers = FindMarkers_UMI(object, ident.1 = opt$type,
                              group.by = opt$ident,
                              assay = "SCT",
                              logfc.threshold = 0.5,
                                 only.pos = T,
                                 test.use = "MAST",
                                 latent.vars = "nFeature_SCT")
    markers$cluster = as.character(opt$type)
    markers$Cell_category = opt$ident
    num = opt$num
    if(args < 10) num = paste0("0",num)
    if(args < 100) num = paste0("0",num)

    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)

    write.csv(markers,paste0(path,arg,"-",opt$ident,"-",num,".",opt$type, ".csv"))
}
