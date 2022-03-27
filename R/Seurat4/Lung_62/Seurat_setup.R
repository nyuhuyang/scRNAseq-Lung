########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.1.1
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
df_samples <- readxl::read_excel("output/20220322/20220322_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()

df_samples = df_samples[!grepl("Invitro", df_samples$condition),]

nrow(df_samples)
#======1.2 load  Seurat =========================
object = readRDS(file = "data/Lung_62_20220322.rds")

table(df_samples$sample.id %in% object$orig.ident)
meta.data = object@meta.data
table(df_samples$sample.id %in% meta.data$orig.ident)
for(i in 1:length(df_samples$sample)){
    cells <- meta.data$orig.ident %in% df_samples$sample[i]
    print(df_samples$sample[i])
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
RowSums = rowSums(object[["RNA"]]@data)
c("SCGB3A2","SCGB3A1","SCGB1A1") %in% rownames(object)[RowSums >= 50]
jpeg(paste0(path,"RowSums.jpeg"), units="in", width=10, height=7,res=600)
hist(log(RowSums+1))
abline(v = log(50))
dev.off()
table(RowSums >= 50)
object = object[RowSums >= 50,]
#======1.6 Performing SCTransform and integration =========================

format(object.size(object),unit = "GB")
options(future.globals.maxSize= object.size(object)*10)
Normal = c("SCT","RNA")[2]
if(Normal == "SCT") object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", 
                        verbose = TRUE)
if(Normal == "RNA"){
    object %<>% NormalizeData()
    object <- FindVariableFeatures(object = object, selection.method = "vst",
                                   num.bin = 20, nfeatures = 2000,
                                   mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
    object %<>% ScaleData(verbose = FALSE)
}


object %<>% RunPCA(verbose = T,npcs = 100)

jpeg(paste0(path,"S1_ElbowPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 100)
dev.off()

saveRDS(object, file = "data/Lung_62_20220322.rds")

#======1.8 UMAP from harmony =========================
DefaultAssay(object) = "SCT"

npcs = 100
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
#system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))

object[["harmony.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                 key = "harmonyUMAP_", assay = DefaultAssay(object))
#object[["harmony.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
#                                                 key = "harmonytSNE_", assay = DefaultAssay(object))

npcs = 100
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
#system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
resolutions = c( 0.01, 0.1, 0.2, 0.5,0.8)
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    Progress(i,length(resolutions))
}

saveRDS(object, file = "data/Lung_62_20220322.rds")


# ======1.8.5 Unimodal UMAP Projection =========================
object = readRDS("data/Lung_62_20220322.rds")
DefaultAssay(object) = "SCT"
object[["raw.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "rawUMAP_", assay = DefaultAssay(object))


object.reference = readRDS("data/Lung_SCT_30_20210831.rds")
object.reference[["orig.umap"]] <- CreateDimReducObject(embeddings = object.reference@reductions[["umap"]]@cell.embeddings,
                                                 key = "origUMAP_", assay = DefaultAssay(object.reference))
object.reference %<>% RunUMAP(reduction = "pca", dims = 1:100,return.model = TRUE)

embedding = object.reference[["orig.umap"]]@cell.embeddings
colnames(embedding) = c("UMAP_1","UMAP_2")
object.reference[["umap"]]@cell.embeddings = embedding
colnames(embedding) =NULL
object.reference[["umap"]]@misc$model$embedding = embedding
#length(setdiff(colnames(object),colnames(object_reference)))
#newCells <- setdiff(colnames(object),colnames(object.reference))
#object.query = object[,newCells]

anchors <- FindTransferAnchors(reference = object.reference, query = object,
                               dims = 1:30, reference.reduction = "pca")
object <- MapQuery(anchorset = anchors,
                     reference = object.reference, query = object,
                     refdata = list(Cell_subtype = "Cell_subtype"),
                     reference.reduction = "pca", reduction.model = "umap")
colnames(object@meta.data) %<>% gsub("predicted.","",.)

length(intersect(colnames(object),colnames(object.reference)))
oldCells <- intersect(colnames(object),colnames(object.reference))
object@meta.data[oldCells,"Cell_subtype"] = object.reference@meta.data[oldCells,"Cell_subtype"]
#object.query = object[,newCells]

meta.data = object.reference@meta.data
meta.data = meta.data[!duplicated(meta.data$Cell_subtype),]
object$Cell_subtype.colors = plyr::mapvalues(object$Cell_subtype,
                                             from = meta.data$Cell_subtype,
                                             to = meta.data$Cell_subtype.colors)
saveRDS(object, file = "data/Lung_62_20220322.rds")


#=======1.9 save SCT only =======================================
format(object.size(object),unit = "GB")

format(object.size(object@assays$RNA),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")
object[['RNA']] <- NULL
object[["SCT"]]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/Lung_SCT_62_20220322.rds")


