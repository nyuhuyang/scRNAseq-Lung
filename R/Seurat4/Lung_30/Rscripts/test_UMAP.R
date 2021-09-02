# test
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","magrittr","sctransform"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

test_df = data.frame(min_dist = rep(2:5/10,each = 5),
                     spread = rep(3:7/5,times = 4))

print(spread <- test_df[args,"spread"])
print(min.dist <- test_df[args,"min_dist"])
file.name = paste0("dist.",min.dist,"_spread.",spread)

object = readRDS(file = "data/Lung_30_20210831.rds")
DefaultAssay(object) = "SCT"
print(length(VariableFeatures(object)))

npcs = 100
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs,min.dist = min.dist,spread = spread)
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
object[[paste0("umap_",file.name)]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                            key = paste0(args,"UMAP_"), assay = DefaultAssay(object))
umap = object@reductions[paste0("umap_",file.name)]
saveRDS(umap, file = paste0(save.path, "/umap_",file.name,".rds"))

meta.data = object@meta.data[,grep("SCT_snn_res",colnames(object@meta.data),invert = TRUE)]
object@meta.data = meta.data
resolutions = c(0.8,seq(1,5, by = 1))

for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    Progress(i,length(resolutions))
}

colnames(object@meta.data) %<>% gsub("SCT_snn_res",file.name,.)

meta.data = object@meta.data[,grep(file.name,colnames(object@meta.data))]
saveRDS(meta.data, file = paste0(save.path, "/meta.data_",file.name,".rds"))
