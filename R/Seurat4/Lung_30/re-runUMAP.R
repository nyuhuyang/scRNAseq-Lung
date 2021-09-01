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


#########################

# ====== load old Lung data=============
Lung = readRDS(file = "data/Lung_SCT_30_20200710.rds")

# ====== load new Lung data=============
object = readRDS(file = "data/Lung_59_20210814.rds")

Lung$barcode = gsub(".*_","",colnames(Lung)) %>% gsub("-1$","",.)

Lung$orig.ident = gsub("-","_",Lung$orig.ident)
Lung$orig.ident %<>% gsub("VU19_D","VU_19_D",.)
table(Lung$orig.ident %in% object$orig.ident)
Lung$barcode = paste0(Lung$orig.ident,"-",Lung$barcode)
Lung %<>% RenameCells(new.names = Lung$barcode)

table(Lung$barcode %in% colnames(object))

meta_data = Lung@meta.data[,c("barcode","cell_types","cell_types.colors","Doublets")]
rownames(meta_data) = meta_data$barcode
meta_data = meta_data[Lung$barcode %in% colnames(object),]
object$barcode = colnames(object)
table(meta_data$barcode %in% object$barcode)
meta.data = object@meta.data
meta.data$cell_types = NULL
meta.data$cell_types.colors = NULL
meta.data %<>% left_join(meta_data,by = "barcode")
meta.data[is.na(meta.data$cell_types),"cell_types"] = "Unknown"

rownames(meta.data) = meta.data$barcode
table(rownames(meta.data) == rownames(object@meta.data))

object@meta.data = meta.data
object %<>% subset(subset = cell_types != "Unknown")
UMAPPlot.1(object,do.print= T, pt.size = 0.5,raster=FALSE, group.by = "cell_types")
remove(Lung);GC()

object[["new.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "newUMAP_", assay = DefaultAssay(object))

object[["new.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
                                             key = "newtSNE_", assay = DefaultAssay(object))

table(colnames(Lung) %in% colnames(object))
table(colnames(object) %in% colnames(Lung))
Lung %<>% subset(cells = colnames(object))

object[["hg19.umap"]] <- CreateDimReducObject(embeddings = Lung@reductions[["umap"]]@cell.embeddings,
                                             key = "hg19UMAP_", assay = DefaultAssay(Lung))

object@reductions$pca = NULL
object@reductions$new.umap = NULL
object@reductions$new.tsne = NULL

#=====================
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
object %<>% FindVariableFeatures(selection.method = "vst",
                                 num.bin = 20,
                                 mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

object %<>% ScaleData
npcs <- 105
object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))

jpeg(paste0(path,"S1_ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = npcs)
dev.off()

object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

meta.data = object@meta.data
meta.data = meta.data[!duplicated(meta.data$cell_types),]
meta.data = meta.data[order(meta.data$cell_types),]

resolutions = c(seq(0.1,5, by = 0.1),seq(2,5, by = 1))
for(i in 1:length(resolutions)){
        object %<>% FindClusters(resolution = resolutions[i])
        Progress(i,length(resolutions))
}

UMAPPlot.1(object,do.print= T, pt.size = 0.5,raster=FALSE,cols=meta.data$cell_types.colors, group.by = "cell_types")
saveRDS(object, file = "data/Lung_30_20210831.rds")

#=======1.9 save SCT only =======================================
format(object.size(object),unit = "GB")

format(object.size(object@assays$RNA),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")
object[['RNA']] <- NULL
object[['integrated']] <- NULL
object[["SCT"]]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/Lung_SCT_30_20210831.rds")

