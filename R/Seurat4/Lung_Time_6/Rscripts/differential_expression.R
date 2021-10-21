library(Seurat)
library(magrittr)
library(tidyr)
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

#==========================
object = readRDS(file = "data/Lung_SCT_time6_20210908.rds")

step = c("resolutions","Analysis 1a","Analysis 1b","Analysis 2")[4]

if(step == "resolutions"){
    opts = paste0("SCT_snn_res.",c(0.01, 0.05, 0.09, 0.1, 0.5, 0.6, 0.7))
    opt = opts[args]
    print(opt)
    
    Idents(object) = opt
    
    markers = FindAllMarkers_UMI(object,
                                 group.by = opt,
                                 assay = "SCT",
                                 logfc.threshold = 0.5,
                                 only.pos = T,
                                 test.use = "MAST",
                                 latent.vars = "nFeature_SCT")
    
    write.csv(markers,paste0(path,args,"_",opt, ".csv"))
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
