########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","magrittr",
                   "tidyr","gplots","MAST"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# Load Seurat
(load(file="data/Lung_24_20190918.Rda"))
Idents(object) = "conditions"
table(Idents(object))
#======1.5 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "conditions")
remove(object);GC()
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
npcs =50
object_list %<>% lapply(function(x) {
        x %<>% RunPCA(features = object.features, verbose = FALSE,npcs = npcs)
})
options(future.globals.maxSize= object.size(object_list)*3)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features, 
                                  verbose = FALSE)
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT", 
                                  anchor.features = object.features,
                                  reference = c(1, 2), reduction = "rpca", 
                                  dims = 1:npcs)
object <- IntegrateData(anchorset = anchors, normalization.method = "SCT",dims = 1:npcs)
save(object, file = paste0("data/Lung_24_20191105.Rda"))
remove(anchors,object_list);GC()
# FindNeighbors
object %<>% RunPCA(npcs = 100, verbose = FALSE)
npcs =100
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
UMAPPlot.1(object = object, group.by="integrated_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T,do.print = T,unique.name = "conditions",
           label.size = 4, repel = T,no.legend = T,
           title = paste("UMAP plot of 24 samples"))

TSNEPlot.1(object, group.by="integrated_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T,do.print = T,unique.name = "conditions",
           label.size = 4, repel = T,no.legend = T,
           title = paste("tSNE plot of 24 samples"))
object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = paste0("data/Lung_24_20191105.Rda"))

DefaultAssay(object) = "integrated"
npcs =100
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
for(i in c(4:16)/10){
        object %<>% FindClusters(resolution = i)
        Idents(object) = paste0("integrated_snn_res.",i)
        TSNEPlot.1(object, group.by=paste0("integrated_snn_res.",i),
                   pt.size = 1,label = T,
                   label.repel = T,do.print = T,
                   label.size = 4, repel = T,no.legend = T,
                   title = paste("resolution =",i))
}

DefaultAssay(object) = "RNA"
npcs =100
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
for(i in c(4:16)/10){
        object %<>% FindClusters(resolution = i)
        Idents(object) = paste0("RNA_snn_res.",i)
        TSNEPlot.1(object, group.by=paste0("RNA_snn_res.",i),
                   pt.size = 1,label = T,
                   label.repel = T,do.print = T,
                   label.size = 4, repel = T,no.legend = T,
                   title = paste("resolution =",i))
}