########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony","magrittr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel(paste0(path,"20210831_scRNAseq_info.xlsx"))
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()

nrow(df_samples)
df_samples$date %<>% gsub(" UTC","",.) %>% as.character()
#======1.2 load  Seurat =========================
object = readRDS(file = "data/Lung_62_20210831.rds")

table(df_samples$sample.id %in% object$orig.ident)
meta.data = object@meta.data
table(df_samples$sample.id %in% meta.data$orig.ident)
for(i in 1:length(df_samples$sample.id)){
    cells <- meta.data$orig.ident %in% df_samples$sample.id[i]
    print(df_samples$sample.id[i])
    print(table(cells))
    meta.data[cells,"regions"] = df_samples$regions[i]
    meta.data[cells,"patient"] = df_samples$patient[i]
    meta.data[cells,"condition"] = as.character(df_samples$condition[i])
    meta.data[cells,"date"] = as.character(df_samples$date[i])
    meta.data[cells,"project"] = as.character(df_samples$project[i])
    meta.data[cells,"Mean.Reads.per.Cell"] = df_samples$mean.reads.per.cell[i]
    meta.data[cells,"Number.of.Reads"] = df_samples$number.of.reads[i]
    meta.data[cells,"Sequencing.Saturation"] = df_samples$sequencing.saturation[i]
}
meta.data$orig.ident %<>% factor(levels = df_samples$sample.id)
table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data
Idents(object) = "orig.ident"
#======1.6 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()

object_list %<>% lapply(function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = object_list)

options(future.globals.maxSize= object.size(object_list)*1.5)
object_list %<>% pblapply(FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = object_list,anchor.features = features,
                                  reference = c(1, 2, 3), reduction = "rpca")

#remove(object_list)
GC()
# this command creates an 'integrated' data assay
object <- IntegrateData(anchorset = anchors)
remove(anchors);GC()
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/Lung_62_20210831.rds")
# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells!

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(object) <- "integrated"

# Run the standard cca workflow for umap & tsne visualization
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = 100, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = 100))
dev.off()
npcs = 100

object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))

object[["cca.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "ccaUMAP_", assay = DefaultAssay(object))
object[["cca.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
                                             key = "ccatSNE_", assay = DefaultAssay(object))

saveRDS(object, file = "data/Lung_62_20210831.rds")


#======1.7 UMAP from raw pca =========================
format(object.size(object),unit = "GB")
DefaultAssay(object) = "RNA"
#object[["SCT]] is too large, around 80 GB
object[['SCT']] = NULL
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)

object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = 100)

jpeg(paste0(path,"S1_ElbowPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 100)
dev.off()

saveRDS(object, file = "data/Lung_62_20210831.rds")

#======1.8 UMAP from harmony =========================
DefaultAssay(object) = "SCT"

npcs = 100
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))

object[["harmony.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                 key = "harmonyUMAP_", assay = DefaultAssay(object))
object[["harmony.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
                                                 key = "harmonytSNE_", assay = DefaultAssay(object))

npcs = 100
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
object %<>% FindClusters(resolution = 0.8)
resolutions = c(seq(0.01,0.09, by = 0.01),seq(0.1,0.9, by = 0.1),seq(1,5, by = 0.1))
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    Progress(i,length(resolutions))
}

saveRDS(object, file = "data/Lung_62_20210831.rds")


#=======1.9 save SCT only =======================================
format(object.size(object),unit = "GB")

format(object.size(object@assays$RNA),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")
object[['RNA']] <- NULL
object[['integrated']] <- NULL
object[["SCT"]]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/Lung_SCT_62_20210831.rds")


