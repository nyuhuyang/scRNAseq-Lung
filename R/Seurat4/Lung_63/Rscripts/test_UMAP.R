# r4.1.1
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

# Need 64GB
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

#=========
test_df = data.frame(min_dist = rep(c(0.1,0.2,0.5),each = 3),
                     spread = rep(c(0.5,1.0,1.5),times = 3),
                     npcs = rep(c(80,90,100,110,120,130),each = 9))
test_df = bind_rows(list(test_df,test_df))
test_df$nfeatures = rep(c(2000,3000),each =54)

print(nfeatures <- test_df[args,"nfeatures"])
print(spread <- test_df[args,"spread"])
print(min.dist <- test_df[args,"min_dist"])
print(npcs <- test_df[args,"npcs"])

file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread)

save.path <- paste0("output/",gsub("-","",Sys.Date()),"/",nfeatures)
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)


object = readRDS(file = "data/Lung_SCT_63_20220408.rds")
object@meta.data = readRDS(file = "data/Lung_63_20220408_meta.data.rds")
DefaultAssay(object) = "SCT"

object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = nfeatures,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
print(length(VariableFeatures(object)))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs,min.dist = min.dist,spread = spread,
                    return.model = TRUE)

UMAPPlot.1(object, group.by = "Cell_subtype",do.print = T,raster=FALSE,no.legend = T,
           title = file.name,file.name = paste0(file.name,".jpeg"),
           save.path = save.path)

object[[paste0("umap_",file.name)]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                            key = paste0(args,"UMAP_"), assay = DefaultAssay(object))
umap = object@reductions[paste0("umap_",file.name)]
saveRDS(umap, file = paste0(save.path, "/umap_",file.name,".rds"))
