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

object = readRDS(file = "data/Lung_SCT_62_20220322.rds")

(step = c("resolutions","cell_types")[1])

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

if(step == "cell_types"){# 32~64GB
    meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")
    table(rownames(object@meta.data) == rownames(meta.data))
    table(object$barcode ==meta.data$barcode)
    object@meta.data = meta.data
    
    df = matrix(c('AT1', 'AT',
                  'AT2', 'AT',
                  'B',   'B',
                  'BC-IC','BC-IC',
                  'C',    'C',
                  'CD4-T','T',
                  'CD8-T','T',
                  'Cr',   'Cr',
                  'DC',   'MPS',
                  'En-a', 'En-a',
                  'En-c', 'En-c',
                  'En-l', 'En-l',
                  'En-SM', 'En-SM',
                  'En-v', 'En-v',
                  'Fb',   'Fb',
                  'G-Muc','G-Muc',
                  'G-Ser','G-Ser',
                  'Gli', 'Gli', 
                  'Ion','Ion',
                  'M', 'MPS',
                  'MC','MC',
                  'ME', 'ME',
                  'Mon','MPS',
                  'NE','NE',
                  'Neu','Neu',
                  'NK','NK',
                  'PC','PC',
                  'Pr','SM+Pr',
                  'S','S',
                  'SCI','SCI',
                  'SM','SM+Pr',
                  'T-NK','T',
                  'T-un','T'),
                ncol = 2, byrow = TRUE,
                dimnames = list(1:33,
                                c("Cell_type", "Cell_group"))) %>%
        as.data.frame()
    DefaultAssay(object) = "SCT"
    object %<>% subset(subset = Cell_subtype != "Un"
                       &  Doublets == "Singlet"
    )
    object$Cell_group = plyr::mapvalues(object$Cell_type,
                                        from = df$Cell_type,
                                        to = df$Cell_group)
    opts = data.frame(ident = c(rep("Cell_subtype",49),
                                rep("Cell_type",31),
                                rep("Cell_group",26), #81~106
                                rep("UMAP_land",20),
                                rep("Family",7),
                                rep("Superfamily",3)),
                      num = c(1:49,
                              1:31,
                              1:26,
                              1:20,
                              1:7,
                              1:3)
                      )
    opt = opts[args,]

    #==========================
    Idents(object) = opt$ident
    opt$type = sort(levels(object))[opt$num]
    print(opt)
    

    markers = FindMarkers_UMI(object, ident.1 = opt$type,
                              group.by = opt$ident,
                              assay = "SCT",
                              min.pct = 0.01,
                              logfc.threshold = 0.1,
                                 only.pos = F,
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

    save_path <- paste0(path,step,"/")
    if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
    
    write.csv(markers,paste0(save_path,arg,"-",opt$ident,"-",num,".",opt$type, ".csv"))
}
