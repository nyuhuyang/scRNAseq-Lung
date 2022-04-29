# r4.1.1
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

# read sample summary list
df_samples <- readxl::read_excel("output/20220408/20220406_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
nrow(df_samples)
#======================================
object = readRDS(file = "data/Lung_SCT_63_20220408.rds")
object@meta.data = readRDS(file = "data/Lung_63_20220408_meta.data.rds")

test_df = data.frame(min_dist = rep(c(0.1,0.2,0.5),each = 3),
                     spread = rep(c(0.5,1.0,1.5),times = 3),
                     npcs = rep(c(80,90,100,110,120,130),each = 9))
test_df = bind_rows(list(test_df,test_df))
test_df$nfeatures = rep(c(3000,2000),each =54)

for(i in c(21,75)){
    nfeatures <- test_df[i,"nfeatures"]
    spread <- test_df[i,"spread"]
    min.dist <- test_df[i,"min_dist"]
    npcs <- test_df[i,"npcs"]
    umap.name = paste0("nfeatures",nfeatures,"npcs",npcs,"dist.",min.dist,"spread.",spread) %>% gsub("\\.","",.)
    print(i);    print(umap.name)
    file.name = paste0("output/20220420/",nfeatures,"/umap_cs",npcs,"_dist.",min.dist,"_spread.",spread,".rds")
    umap  = readRDS(file.name)[[1]]
    umap@key = paste0(umap.name,"_")
    colnames(umap@cell.embeddings) = paste0(umap.name,"_",1:2)
    object[[paste0("umap_",umap.name)]] <- umap
    Progress(i,108)
}

SCT_snn <- grep("SCT_snn_res",colnames(object@meta.data), value = T)
for(cl in SCT_snn) object@meta.data[,cl] %<>% as.character() %<>% as.integer()

