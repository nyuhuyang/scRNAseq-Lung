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

object = readRDS(file = "data/Lung_SCT_30_20210831.rds")

(step = c("resolutions","cell_types",
         "DEG analysis 1","DEG analysis 2",
         "DEG analysis 3-option 1",
         "DEG analysis 3-option 2","ASE",
         "TASCs lvl")[3])

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


if(step == "DEG analysis 1"){
    # DEG analysis 1 - between TASC and related cell populations (with volcanos; cells in the same group from different samples mixed):
    meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")
    table(rownames(object@meta.data) == rownames(meta.data))
    table(object$barcode == meta.data$barcode)
    object@meta.data = meta.data
    
    DefaultAssay(object) = "SCT"
    Region = switch(as.character(args < 11),
                     "TRUE" ="distal",
                     "FALSE"= c("distal","terminal"))
    Regions_name = ifelse(args < 11, "_D", "_DT")

    Idents_list = list(ident1 = list("TASC",
                                     "TASC",
                                     "TASC",
                                     "TASC",
                                     "TASC",
                                     "TASC",
                                     "S-Muc",
                                     "S-Muc",
                                     "S-Muc",
                                     "S-Muc",
                                     c("S1", "S-Muc","TASC")),
                       ident2 = list("S1",
                                     c("S1", "S-Muc"),
                                     "G-Ser",
                                     c("G-Ser","G-Muc"),
                                     "AT2",
                                     c("AT1","AT2"),
                                     "S1",
                                     c("TASC","S1"),
                                     "G-Muc",
                                     c("G-Muc","G-Ser"),
                                     c("G-Muc","G-Ser")))
    
    (ident1 <- Idents_list$ident1[[ args %% 11 +1]])# 0:21
    (ident2 <- Idents_list$ident2[[args %% 11 +1]])# 0:21
    
    object %<>% subset(subset = Cell_subtype %in% c(ident1,ident2)
                       &  Doublets == "Singlet"
                       & Regions %in% Region
    )
    
    #==========================
    Idents(object) = "Cell_subtype"
    markers = FindMarkers_UMI(object, ident.1 = ident1,
                              ident.2 = ident2,
                              group.by = "Cell_subtype",
                              assay = "SCT",
                              logfc.threshold = 0.1,
                              min.pct = 0.01,
                              only.pos = F,
                              test.use = "MAST",
                              latent.vars = "nFeature_SCT")
    markers$cluster = paste(paste(ident1,collapse = "+"),
                            "vs.",
                            paste(ident2,collapse = "+"))
    markers$gene = rownames(markers)

    arg = args
    if(args < 10) arg = paste0("0",arg)

    save_path <- paste0(path,step,"/")
    if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
    
    file.name = paste0(paste(ident1,collapse = "_"),
                      "-vs-",
                      paste(ident2,collapse = "_"))
    write.csv(markers,paste0(save_path,arg,"-",file.name,Regions_name, ".csv"))
}


if(step == "DEG analysis 2"){
    # DEG analysis 2 – between different subsets of TASCs
    meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")
    table(rownames(object@meta.data) == rownames(meta.data))
    table(object$barcode == meta.data$barcode)
    object@meta.data = meta.data
    
    DefaultAssay(object) = "SCT"
    
    
    Idents_list = list(ident1 = list("SCGB1A1-lo TASC",
                                     "SCGB3A2+ TASC",
                                     "SCGB1A1-hi TASC",
                                     "SCGB1A1-hi TASC",
                                     "SCGB1A1-lo TASC",
                                     "SCGB1A1-lo TASC",
                                     "SCGB3A2- TASC",
                                     "SCGB3A2- TASC",
                                     "SCGB3A2+ TASC",
                                     "SCGB3A2+ TASC"),
                       ident2 = list("SCGB1A1-hi TASC",
                                     "SCGB3A2- TASC",
                                     "S1",
                                     c("S1", "S-Muc"),
                                     "AT2",
                                     c("AT1","AT2"),
                                     "S1",
                                     c("S1","S-Muc"),
                                     "AT2",
                                     c("AT1","AT2")))
    group.by.list = c("Cell_subtype.SCGB1A1",
                      "Cell_subtype.SCGB3A2",
                      rep("Cell_subtype.SCGB1A1",4),
                      rep("Cell_subtype.SCGB3A2",4))
    (ident1 <- Idents_list$ident1[[args]])
    (ident2 <- Idents_list$ident2[[args]])
    c(group.by <- group.by.list[args])
    
    object %<>% subset(subset = Doublets == "Singlet")

    #==========================
    Idents(object) = group.by
    markers = FindMarkers_UMI(object, ident.1 = ident1,
                              ident.2 = ident2,
                              group.by = group.by,
                              assay = "SCT",
                              logfc.threshold = 0.1,
                              min.pct = 0.01,
                              only.pos = F,
                              test.use = "MAST",
                              latent.vars = "nFeature_SCT")
    markers$cluster = paste(ident1,
                            "vs.",
                            paste(ident2,collapse = "+"))
    markers$gene = rownames(markers)
    
    arg = args
    if(args < 10) arg = paste0("0",arg)
    
    save_path <- paste0(path,step,"/")
    if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
    
    file.name = paste0(paste(ident1,collapse = "_"),
                       "-vs-",
                       paste(ident2,collapse = "_"))
    write.csv(markers,paste0(save_path,arg,"-",file.name, ".csv"))
}


if(step == "DEG analysis 3-option 1"){
    meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")
    table(rownames(object@meta.data) == rownames(meta.data))
    table(object$barcode ==meta.data$barcode)
    object@meta.data = meta.data
    
    DefaultAssay(object) = "SCT"
    object %<>% subset(subset = Cell_subtype != "Un"
                       &  Doublets == "Singlet"
    )
    Cell_types <- c("Cell_subtype","Cell_type","UMAP_land","Family","Superfamily")
    categories = lapply(Cell_types, function(Cell_type){
        unique(object@meta.data[,Cell_type])
    })
    names(categories) = Cell_types
    category_df = data.frame("category" = c(rep(names(categories[1]),49),
                                rep(names(categories[2]),33),
                                rep(names(categories[3]),20),
                                rep(names(categories[4]),7),
                                rep(names(categories[5]),3)),
                      "type" = c(categories[[1]],
                              categories[[2]],
                              categories[[3]],
                              categories[[4]],
                              categories[[5]])
                      )
    category_df = category_df[!duplicated(category_df$type),]
    category_df = category_df[category_df$type != "Un",]
    
    Idents_list = list(ident1 = list("distal",
                                     "proximal",
                                     "proximal",
                                     c("distal", "terminal"),
                                     "distal",
                                     c("proximal", "terminal"),
                                     "terminal",
                                     "distal",
                                     "terminal",
                                     c("proximal", "distal"),
                                     "COPD",
                                     "distal"),
                       ident2 = list("proximal",
                                     "distal",
                                     c("distal", "terminal"),
                                     "proximal",
                                     c("proximal","terminal"),
                                     "distal",
                                     "distal",
                                     "terminal",
                                     c("proximal","distal"),
                                     "terminal",
                                     "distal",
                                     "COPD"))
    i = ceiling((args/75) %% 12) # 75 cell types, 12 pairs
    if(args == 900) i = 12
    print(ident1 <- Idents_list$ident1[[i]])
    print(ident2 <- Idents_list$ident2[[i]])
    
    k = ((args-1) %% 75)+1
    print(category <- category_df$category[k])
    print(type <- as.character(category_df$type[k]))
    
    #==========================
    Idents(object) = category
    object %<>% subset(idents = type)
    
    Idents(object) = "Regions"
    
    markers = FindMarkers_UMI(object, ident.1 = ident1,
                              ident.2 = ident2,
                              group.by = "Regions",
                              assay = "SCT",
                              min.pct = 0.01,
                              logfc.threshold = 0.1,
                              only.pos = T,
                              test.use = "wilcox")
    change_name <- function(ident) {
        switch (ident,
        "distal" = "D",
        "terminal" = "T",
        "proximal" = "P",
        "COPD" = "COPD",
        "distal+terminal" = "D+T",
        "proximal+terminal" = "P+T",
        "proximal+distal" = "P+D"
        )
    }
    ident1 = change_name(paste(ident1,collapse = "+"))
    ident2 = change_name(paste(ident2,collapse = "+"))
    
    markers$DE_pairs = paste(ident1, "vs", ident2)
    markers$type = type
    markers$Cell_category = category
    
    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)
    if(args < 1000) arg = paste0("0",arg)
    
    save_path <- paste0(path,step,"/")
    if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
    
    write.csv(markers,paste0(save_path,arg,"-",ident1,"-vs-",ident2,".",type, ".csv"))
}


if(step == "DEG analysis 3-option 2"){
    meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")
    table(rownames(object@meta.data) == rownames(meta.data))
    table(object$barcode ==meta.data$barcode)
    object@meta.data = meta.data
    
    DefaultAssay(object) = "SCT"
    object %<>% subset(subset = Cell_subtype != "Un"
                       &  Doublets == "Singlet"
    )
    Cell_types <- c("Cell_subtype","Cell_type","UMAP_land","Family","Superfamily")
    categories = lapply(Cell_types, function(Cell_type){
        unique(object@meta.data[,Cell_type])
    })
    names(categories) = Cell_types
    category_df = data.frame("category" = c(rep(names(categories[1]),49),
                                            rep(names(categories[2]),33),
                                            rep(names(categories[3]),20),
                                            rep(names(categories[4]),7),
                                            rep(names(categories[5]),3)),
                             "type" = c(categories[[1]],
                                        categories[[2]],
                                        categories[[3]],
                                        categories[[4]],
                                        categories[[5]])
    )
    category_df = category_df[!duplicated(category_df$type),]
    category_df = category_df[category_df$type != "Un",]
    
    meta.data = meta.data[!duplicated(meta.data$orig.ident),]
    meta.data$orig.ident %<>% as.character()
    meta.data$Regions %<>% as.character()
    
    all_samples = as.character(meta.data$orig.ident)
    
    Idents_list = list(ident1 = list("distal",
                                     "proximal",
                                     "proximal",
                                     c("distal", "terminal"),
                                     "distal",
                                     c("proximal", "terminal"),
                                     "terminal",
                                     "distal",
                                     "terminal",
                                     c("proximal", "distal"),
                                     "COPD",
                                     "distal"),
                       ident2 = list("proximal",
                                     "distal",
                                     c("distal", "terminal"),
                                     "proximal",
                                     c("proximal","terminal"),
                                     "distal",
                                     "distal",
                                     "terminal",
                                     c("proximal","distal"),
                                     "terminal",
                                     "distal",
                                     "COPD"))
    i = ceiling((args/75) %% 12) # 75 cell types, 12 pairs
    if(args == 900) i = 12
    print(ident1 <- Idents_list$ident1[[i]])
    print(ident2 <- Idents_list$ident2[[i]])
    
    k = ((args-1) %% 75)+1
    print(category <- category_df$category[k])
    print(type <- as.character(category_df$type[k]))
    
    #==========================
    Idents(object) = category
    object = subset(object, idents = type)
    
    Idents(object) = "orig.ident"
    exp = AverageExpression(object,assays = "SCT",verbose = T) %>% .$SCT
    exp = log1p(exp)
    exp %<>% as.data.frame()
    # fill up missing data with 0
    if(ncol(exp) < 30) {
        missing_sample = all_samples[!all_samples %in% colnames(exp)]
        exp[,missing_sample] = 0
    }
    object = CreateSeuratObject(exp,min.cells = 0,names.delim = "-",min.features = 0)
    object$Regions = plyr::mapvalues(colnames(object),
                                         from = meta.data$orig.ident,
                                         to = meta.data$Regions)
    markers = FindMarkers_UMI(object, ident.1 = ident1,
                              ident.2 = ident2,
                              group.by = "Regions",
                              assay = "RNA",
                              min.pct = 0,
                              logfc.threshold = 0,
                              only.pos = T,
                              test.use = "wilcox")
    change_name <- function(ident) {
        switch (ident,
                "distal" = "D",
                "terminal" = "T",
                "proximal" = "P",
                "COPD" = "COPD",
                "distal+terminal" = "D+T",
                "proximal+terminal" = "P+T",
                "proximal+distal" = "P+D"
        )
    }
    ident1 = change_name(paste(ident1,collapse = "+"))
    ident2 = change_name(paste(ident2,collapse = "+"))
    
    markers$DE_pairs = paste(ident1, "vs", ident2)
    markers$type = type
    markers$Cell_category = category
    
    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)

    save_path <- paste0(path,step,"/")
    if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
    
    write.csv(markers,paste0(save_path,arg,"-",ident1,"-vs-",ident2,".",type, ".csv"))
}


if(step == "ASE"){
    meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")
    table(rownames(object@meta.data) == rownames(meta.data))
    table(object$barcode ==meta.data$barcode)
    object@meta.data = meta.data
    
    Region = switch(as.character(args < 11),
                     "TRUE" ="distal",
                     "FALSE"= c("distal","terminal"))
    Regions_name = ifelse(args < 11, "D", "DT")
    DefaultAssay(object) = "SCT"
    object %<>% subset(subset = Cell_subtype != "Un"
                       & Family == "ASE"
                       & Regions %in% Region
                       &  Doublets == "Singlet"
    )
    Cell_types <- c("BC","C-s","C1","H","IC","Ion","NE","p-C","S-Muc","S1","TASC")
    (cell_type = Cell_types[args %% 11 +1]) # 0:21
    #==========================
    Idents(object) = "Cell_subtype"

    markers = FindMarkers_UMI(object, ident.1 = cell_type,
                              group.by = "Cell_subtype",
                              assay = "SCT",
                              min.pct = 0.01,
                              logfc.threshold = 0.1,
                              only.pos = F,
                              test.use = "MAST",
                              latent.vars = "nFeature_SCT")

    arg = args
    if(args < 10) arg = paste0("0",arg)

    save_path <- paste0(path,step,"/")
    if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
    
    write.csv(markers,paste0(save_path,arg,"-",cell_type,"_vs_ASE_",Regions_name, ".csv"))
}


if(step == "TASCs lvl"){
    opts = data.frame(gene = c(rep("SCGB1A1",8),
                                rep("SCGB3A2",3)),
                      lvl = c(0:7,
                              2,1,0)
    )
    opt = opts[args,]
    
    meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")
    table(rownames(object@meta.data) == rownames(meta.data))
    table(object$barcode ==meta.data$barcode)
    object@meta.data = meta.data
    
    DefaultAssay(object) = "SCT"
    TASC <- subset(object, subset = Cell_subtype %in% "TASC" & 
                       Doublets %in% "Singlet" &
                       Regions %in% c("distal","terminal"))
    
    TASC_exp = FetchData(TASC, vars = opt$gene)
    TASC_exp[,opt$gene] = plyr::mapvalues(as.character(TASC_exp[,opt$gene] > opt$lvl),
                                       from = c("TRUE","FALSE"),
                                       to = paste(opt$gene,paste(c("high","low"),opt$lvl)))
    print(table(TASC_exp[,opt$gene] ))
    TASC@meta.data %<>% cbind(TASC_exp)
    #==========================
    Idents(TASC) = opt$gene
    
    markers = FindMarkers_UMI(TASC, ident.1 = paste(opt$gene,"high",opt$lvl),
                              group.by = opt$gene,
                              assay = "SCT",
                              min.pct = 0.01,
                              logfc.threshold = 0.05,
                              only.pos = F,
                              test.use = "MAST",
                              latent.vars = "nFeature_SCT")
    markers$cluster = paste(opt$gene,"high vs low at lvl",opt$lvl)
    markers$gene = rownames(markers)
    write.csv(markers,paste0(path,args,"-",paste(opt$gene,"high vs low at",opt$lvl),"_TACS.csv"))
}
