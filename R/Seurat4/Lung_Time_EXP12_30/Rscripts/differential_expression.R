library(Seurat)
library(magrittr)
library(tidyr)
library(dplyr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/EXP12_30/")
if(!dir.exists(path)) dir.create(path, recursive = T)


set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

#==========================
object = readRDS(file = "data/Lung_time15_SoupX_20230129.rds")
meta.data <- readRDS("output/Lung_time15_metadata_20220523_v2.rds")
if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data <- meta.data
}
step = c("resolutions","DGEs","DGEs_group","DGEs_group2","Analysis 1a","Analysis 1b","Analysis 2",
         "DGEs_pseudobulk")[8]

if(step == "resolutions"){ # Need 64GB 
    opts = data.frame(ident = c(rep("SCT_snn_res.0.3",21),
                                rep("SCT_snn_res.0.5",26),
                                rep("SCT_snn_res.0.8",33),
                                rep("SCT_snn_res.1",35),
                                rep("SCT_snn_res.2",49),
                                rep("SCT_snn_res.5",92)),
                      num = c(0:20,
                              0:25,
                              0:32,
                              0:34,
                              0:48,
                              0:91)
                      )

    opt = opts[args,] #1-202
    print(opt)
    #==========================
    Idents(object) = opt$ident

    markers = FindMarkers_UMI(object, ident.1 = opt$num,
                              group.by = opt$ident,
                              assay = "SCT",
                              min.pct = 0.1,
                              logfc.threshold = 0.25,
                                 only.pos = TRUE,
                                 test.use = "MAST",
                                 latent.vars = "nFeature_SCT"
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

if(step == "DGEs"){ # Need 64GB 
    meta.data[!duplicated(meta.data$`cell state`),
              c("cell state","cell type","cell group","Stage","Path")] %>% 
        tidyr::pivot_longer(everything(),names_to = "ident") %>% 
        distinct(ident,value) %>%
        arrange(ident,value) -> opts
     
    opt = opts[args,] # 1-83
    
    print(opt)
    #==========================
    Idents(object) = opt$ident
    
    markers = FindMarkers_UMI(object, ident.1 = opt$value,
                              group.by = opt$ident,
                              assay = "SCT",
                              min.pct = 0.1,
                              logfc.threshold = 0.1,
                              only.pos = TRUE,
                              test.use = "MAST",
                              latent.vars = "nFeature_SCT"
    )
    markers$category = opt$ident
    markers$cluster = opt$value
    
    arg = args
    if(args < 10) arg = paste0("0",arg)

    write.csv(markers,paste0(path,arg,"_",opt$ident,"_",opt$value, ".csv"))
}

if(step == "DGEs_group"){ # Need 64GB 
    meta.data[!duplicated(meta.data$`cell state`),
              c("cell state","cell type","cell group","Stage","Path")] %>% 
        rename_with(~ gsub(" ", "_", .x, fixed = TRUE)) %>%
        tidyr::pivot_longer(!cell_state,
                            names_to = "ident") %>% 
        arrange(cell_state) -> opts
    
    opt = opts[args,] # 1-220
    
    print(opt)
    #==========================
    colnames(object@meta.data) %<>% gsub(" ", "_", ., fixed = TRUE)
    object <- object[,object@meta.data[, opt$ident] == opt$value]
    Idents(object) = "cell state"
    
    markers = FindMarkers_UMI(object, ident.1 = opt$cell_state,
                              group.by = "cell_state",
                              assay = "SCT",
                              min.pct = 0.1,
                              logfc.threshold = 0.1,
                              only.pos = TRUE,
                              test.use = "MAST",
                              latent.vars = "nFeature_SCT"
    )
    markers$category = opt$ident
    markers$cluster = opt$value
    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)
    
    path <- paste0(path,step,"/")
    if(!dir.exists(path)) dir.create(path, recursive = T)
    write.csv(markers,paste0(path,arg,"_",opt$cell_state,"_in_",opt$ident,"_",opt$value, ".csv"))
}

if(step == "DGEs_group2"){ # Need 128GB 
    meta.data[!duplicated(meta.data$Stage),c("Stage","Path")] %>%
        tidyr::pivot_longer(!Stage,
                            names_to = "ident") %>% 
        arrange(Stage) -> opts
    
    opt = opts[args,] # 1-8
    
    print(opt)
    #==========================
    object %<>% subset(subset = Path == opt$value)
    Idents(object) = "Stage"
    
    markers = FindMarkers_UMI(object, ident.1 = opt$Stage,
                              group.by = "Stage",
                              assay = "SCT",
                              min.pct = 0.1,
                              logfc.threshold = 0.1,
                              only.pos = TRUE,
                              test.use = "MAST",
                              latent.vars = "nFeature_SCT"
    )
    markers$gene = rownames(markers)
    markers$Stage = opt$Stage
    markers$category = opt$ident
    markers$cluster = opt$value
    arg = args
    if(args < 10) arg = paste0("0",arg)
    
    path <- paste0(path,step,"/")
    if(!dir.exists(path)) dir.create(path, recursive = T)
    write.csv(markers,paste0(path,arg,"_",opt$Stage,"_in_",opt$ident,"_",opt$value, ".csv"))
}

if(step == "Analysis 1a"){# 16~32GB
    DefaultAssay(object) = "SCT"
    object %<>% subset(Doublets == "Singlet")
    
    opts = data.frame(ident = c(rep("Cell_subtype",13),
                                rep("Cell_type",7)),
                      type = c(sort(unique(object$Cell_subtype)),
                              sort(unique(object$Cell_type)))
    )
    opts = opts[!duplicated(opts$type),]
    opt = opts[args,]
    
    #==========================
    Idents(object) = opt$ident
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

    arg = args
    if(args < 10) arg = paste0("0",arg)

    save_path <- paste0(path,step,"/")
    if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
    
    write.csv(markers,paste0(save_path,arg,"-",opt$ident,".",opt$type, ".csv"))
}


if(step == "Analysis 1b"){# 16~32GB
    DefaultAssay(object) = "SCT"
    object %<>% subset(subset = Cell_subtype != "Un"
                       &  Doublets == "Singlet"
    )
    object$day %<>% factor(levels = paste0("D",c(0,3,7,14,21,28)))
    df <- as.data.frame.matrix(table(object$Cell_subtype, object$day))
    df %<>% tibble::rownames_to_column(var = "type") %>% 
        pivot_longer(cols = starts_with("D"),
                         names_to = "day",
                         values_to = "cell_number")
    df$ident = "Cell_subtype"
    df1 <- as.data.frame.matrix(table(object$Cell_type, object$day))
    df1 %<>% tibble::rownames_to_column(var = "type") %>% 
        pivot_longer(cols = starts_with("D"),
                     names_to = "day",
                     values_to = "cell_number")
    df1$ident = "Cell_type"
    opts = rbind(df, df1)
    opts$type_day = paste0(opts$type,"_",opts$day)
    opts = subset(opts, cell_number >= 3) 
    opts = opts[!duplicated(opts$type_day),]
    opt = opts[args,]
    
    #==========================
    object$type_day = paste0(object@meta.data[,opt$ident],"_",object$day)
    object %<>% subset(day %in% opt$day)
    
    Idents(object) = "type_day"
    print(opt)
    
    markers = FindMarkers_UMI(object, ident.1 = opt$type_day,
                              group.by = "type_day",
                              assay = "SCT",
                              logfc.threshold = 0.1,
                              only.pos = T,
                              test.use = "MAST",
                              latent.vars = "nFeature_SCT")
    markers$cluster = as.character(opt$type_day)
    markers$Cell_category = opt$ident
    
    arg = args
    if(args < 10) arg = paste0("0",arg)
    
    save_path <- paste0(path,step,"/")
    if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
    
    write.csv(markers,paste0(save_path,arg,"-",opt$ident,".",opt$type, ".csv"))
}



if(step == "Analysis 2"){# 16~32GB
        DefaultAssay(object) = "SCT"
        object %<>% subset(Doublets == "Singlet")
        
        category_df = data.frame(category = c(rep("Cell_subtype",13),
                                    rep("Cell_type",7)),
                          type = c(sort(unique(object$Cell_subtype)),
                                   sort(unique(object$Cell_type)))
        )
        category_df %<>% rbind(c("Cell_type","BC"))
        rownames(category_df) = NULL
        Idents_list = list(ident1 = list("D3",
                                         "D7",
                                         "D7",
                                         "D7",
                                         "D14",
                                         "D14",
                                         "D14",
                                         "D21",
                                         "D21",
                                         "D21",
                                         "D28",
                                         "D28",
                                         "D28"),
                           ident2 = list("D0",
                                         "D3",
                                         "D0",
                                         c("D0","D3"),
                                         "D0",
                                         "D7",
                                         c("D0","D3","D7"),
                                         "D0",
                                         "D14",
                                         c("D0", "D3","D7","D14"),
                                         "D0",
                                         "D21",
                                         c("D0", "D3","D7","D14","D21")))
        
        i = ceiling((args/21) %% 13) # 21 cell types, 13 pairs
        if(args == 273) i = 13
        print(ident1 <- Idents_list$ident1[[i]])
        print(ident2 <- Idents_list$ident2[[i]])
        
        k = ((args-1) %% 21)+1
        print(category <- category_df$category[k])
        print(type <- as.character(category_df$type[k]))
        #==========================
        Idents(object) = category
        object %<>% subset(idents = type)
        
        Idents(object) = "day"
        
        ident1 = ident1[ident1 %in% object$day]
        ident2 = ident2[ident2 %in% object$day]
        
        markers = FindMarkers_UMI(object, ident.1 = ident1,
                                  ident.2 = ident2,
                                  group.by = "day",
                                  assay = "SCT",
                                  min.pct = 0,
                                  logfc.threshold = 0.1,
                                  test.use = "MAST",
                                  latent.vars = "nFeature_SCT")
        markers$category = category
        markers$type = type
        ident1 = paste(ident1,collapse = "+")
        ident2 = paste(ident2,collapse = "+")
        markers$day_pairs = paste(ident1, "vs", ident2)
        
        arg = args
        if(args < 10) arg = paste0("0",arg)
        if(args < 100) arg = paste0("0",arg)
        
        save_path <- paste0(path,step,"/")
        if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
        
        write.csv(markers,paste0(save_path,arg,"-",ident1,".",type, ".csv"))
}

if(step == "DGEs"){ # Need 64GB 
    meta.data[!duplicated(meta.data$`cell state`),
              c("cell state","cell type","cell group","Stage","Path")] %>% 
        tidyr::pivot_longer(everything(),names_to = "ident") %>% 
        distinct(ident,value) %>%
        arrange(ident,value) -> opts
    
    opt = opts[args,] # 1-83
    
    print(opt)
    #==========================
    Idents(object) = opt$ident
    
    markers = FindMarkers_UMI(object, ident.1 = opt$value,
                              group.by = opt$ident,
                              assay = "SCT",
                              min.pct = 0.1,
                              logfc.threshold = 0.1,
                              only.pos = TRUE,
                              test.use = "MAST",
                              latent.vars = "nFeature_SCT"
    )
    markers$category = opt$ident
    markers$cluster = opt$value
    
    arg = args
    if(args < 10) arg = paste0("0",arg)
    
    write.csv(markers,paste0(path,arg,"_",opt$ident,"_",opt$value, ".csv"))
}

if(step == "DGEs_pseudobulk"){ # Need 16GB 
    source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Libra.R")
    
    meta.data[!duplicated(meta.data$`cell state`),
              c("cell state","cell type","cell group","Stage","Path")] %>% 
        tidyr::pivot_longer(everything(),names_to = "ident") %>% 
        distinct(ident,value) %>%
        arrange(ident,value) -> opts
    
    opt = opts[args,] # 1-83
    
    print(opt)
    #==========================
    Idents(object) = opt$ident
    object$cell_type <- "All_celles"
    
    markers <- run_de.1(object, de_family = 'pseudobulk', cell_type_col = "cell_type", 
                        label_col = opt$ident,replicate_col = "orig.ident",
                        ident.1 = opt$value,
                        min_cells = 3,min_reps = 1,
                        de_method = 'edgeR', de_type = 'LRT')

    markers$category = opt$ident

    arg = args
    if(args < 10) arg = paste0("0",arg)
    
    write.csv(markers,paste0(path,arg,"_",opt$ident,"_",opt$value, ".csv"))
}
