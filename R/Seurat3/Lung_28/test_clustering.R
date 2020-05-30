invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","magrittr"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#==== load Seurat ============
(load(file = "data/Lung_28_20200103.Rda"))
object[["cell_types.colors"]] = object[["cell.types.colors"]]
#==== test RNA slot ===========
rm_COPD = T
if(rm_COPD == T) {
        Idents(object) = "conditions"
        object %<>% subset(idents = "COPD",invert = T)
}
DefaultAssay(object) = "RNA"
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
head(VariableFeatures(object), 20)
length(VariableFeatures(object))
object %<>% ScaleData(features = VariableFeatures(object))
npcs =100
object@reductions$pca = NULL
object@reductions$tsne = NULL
object@reductions$umap = NULL
object %<>% RunPCA(features = VariableFeatures(object),verbose =F,npcs = 100)
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
#==== test SCT slot 3000/2000 ===========
test_2000 = T
rm_COPD = T
if(rm_COPD == T) {
        Idents(object) = "conditions"
        object %<>% subset(idents = "COPD",invert = T)
        table(Idents(object))
}
DefaultAssay(object) = "SCT"
if(test_2000 == T){ 
        object <- FindVariableFeatures(object = object, selection.method = "vst",
                                       num.bin = 20,
                                       mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
}
head(VariableFeatures(object), 20)
length(VariableFeatures(object))
if(test_2000 == T) object %<>% ScaleData(features = VariableFeatures(object))
object@reductions$pca = NULL
object@reductions$tsne = NULL
object@reductions$umap = NULL
npcs = 100
object %<>% RunPCA(features = VariableFeatures(object),verbose =F,npcs = npcs)
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

#============ save cell_types plots =========
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
object$cell_types <- plyr::mapvalues(object$cell.types,
        from = df_cell_types$`Cell types`,
        to = df_cell_types$Abbreviation)
Idents(object) = "cell_types"
object %<>% sortIdent()
object %<>% AddMetaColor(label= "cell_types", colors = c(Singler.colors,Singler.colors))

TSNEPlot.1(object, group.by="cell_types",
           cols = ExtractMetaColor(object),
           pt.size = 1,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = T,label.size = 4, title = paste(length(unique(object$orig.ident)),
                                                       "samples No Integration on",
                                                       DefaultAssay(object)),
           do.print = T,do.return = F)
UMAPPlot.1(object, group.by="cell_types",
           cols = ExtractMetaColor(object),
           pt.size = 1,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = T,label.size = 4, title = paste(length(unique(object$orig.ident)),
                                                       "samples No Integration on",
                                                       DefaultAssay(object)),
           do.print = T,do.return = F)
file.rename(paste0(path,"TSNEPlot_object_cell_types.jpeg"), 
            paste0(path,"TSNEPlot_object_",length(unique(object$orig.ident)),
                   "_",DefaultAssay(object),
                   "_",length(VariableFeatures(object)),"_cell_types_Label.jpeg"))
file.rename(paste0(path,"UMAPPlot_object_cell_types.jpeg"), 
            paste0(path,"UMAPPlot_object_",length(unique(object$orig.ident)),
                   "_",DefaultAssay(object),
                   "_",length(VariableFeatures(object)),"_cell_types_Label.jpeg"))

#============ save SCINA plots =========
(load(file = paste0("output/SCINA_20200110.Rda")))
names(results$cell_labels) = colnames(results$probabilities)
object@meta.data$SCINA = results$cell_labels[rownames(object@meta.data)]

df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
object$SCINA %<>% plyr::mapvalues(
        from = df_cell_types$`Cell types`,
        to = df_cell_types$Abbreviation)
Idents(object) = "SCINA"
object %<>% AddMetaColor(label= "SCINA", colors = c(Singler.colors,Singler.colors))

TSNEPlot.1(object, group.by="SCINA",
           cols = ExtractMetaColor(object),
           pt.size = 1,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = T,label.size = 4, title = paste(length(unique(object$orig.ident)),
                                                       "samples No Integration on",
                                                       DefaultAssay(object)),
           do.print = T,do.return = F)
UMAPPlot.1(object, group.by="SCINA",
           cols = ExtractMetaColor(object),
           pt.size = 1,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = T,label.size = 4, title = paste(length(unique(object$orig.ident)),
                                                       "samples No Integration on",
                                                       DefaultAssay(object)),
           do.print = T,do.return = F)
file.rename(paste0(path,"TSNEPlot_object_SCINA.jpeg"), 
            paste0(path,"TSNEPlot_object_",length(unique(object$orig.ident)),
                   "_",DefaultAssay(object),
                   "_",length(VariableFeatures(object)),"_SCINA.jpeg"))
file.rename(paste0(path,"UMAPPlot_object_SCINA.jpeg"), 
            paste0(path,"UMAPPlot_object_",length(unique(object$orig.ident)),
                   "_",DefaultAssay(object),
                   "_",length(VariableFeatures(object)),"_SCINA.jpeg"))
#========
saveRDS(object@reductions, file = paste0(path,"reductions_",length(unique(object$conditions)),
                                         "_",DefaultAssay(object),
                                         "_",length(VariableFeatures(object)),".rds"))
object@reductions = readRDS(file = paste0(path,"reductions_",length(unique(object$conditions)),
                                          "_",DefaultAssay(object),
                                          "_",length(VariableFeatures(object)),".rds"))
readRDS(file = "output/20200114/reductions_4_SCT_2000.rds")

#==== test ICA slot ===========
(load(file = "data/Lung_28_20200116.Rda"))

Idents(object) = "cell_types"
DefaultAssay(object) = "SCT"

head(VariableFeatures(object), 20)
length(VariableFeatures(object))
object@reductions$pca = NULL
object@reductions$tsne = NULL
object@reductions$umap = NULL
npcs = 100
object %<>% RunICA(features = VariableFeatures(object),verbose =T,nics = npcs)

jpeg(paste0(path,"ICAPlot_cell_types_2_3.jpeg"), units="in", width=10, height=7,res=600)
DimPlot(object,dims = c(2,3), reduction = 'ica', cols = ExtractMetaColor(object),
        pt.size = 1,label = T) + NoLegend()
dev.off()

a <- seq(1,97, by = 6)
b <- a+5
a <- 95
for(i in seq_along(a)){
        jpeg(paste0(path,"DimHeatmap_ica_",a[i],"_",b[i],".jpeg"), units="in", width=10, height=7,res=600)
        DimHeatmap(object, dims = a[i]:b[i],nfeatures = 30,reduction = "ica")
        dev.off() 
}
npcs = 100
object %<>% FindNeighbors(reduction = "ica",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunUMAP(reduction = "ica", dims = 1:npcs)
object %<>% RunTSNE(reduction = "ica", dims = 1:npcs)


df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
object$cell_types <- plyr::mapvalues(object$cell.types,
                                     from = df_cell_types$`Cell types`,
                                     to = df_cell_types$Abbreviation)
Idents(object) = "cell_types"
object %<>% sortIdent()
object %<>% AddMetaColor(label= "cell_types", colors = c(Singler.colors,Singler.colors))


TSNEPlot.1(object, group.by="cell_types",
           cols = ExtractMetaColor(object),
           pt.size = 0.5,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = T,label.size = 3, 
           title = "tSNE Plot based on ICA for 28 samples",
           do.print = T,do.return = F)
UMAPPlot.1(object, group.by="cell_types",
           cols = ExtractMetaColor(object),
           pt.size = 0.5,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = T,label.size = 3, 
           title = "UMAP Plot based on ICA for 28 samples",
           do.print = T,do.return = F)

TSNEPlot.1(object, group.by="SCINA",
           cols = ExtractMetaColor(object),
           pt.size = 0.5,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = T,label.size = 3, 
           title = "tSNE Plot based on ICA for 28 samples SCINA",
           do.print = T,do.return = F)
UMAPPlot.1(object, group.by="SCINA",
           cols = ExtractMetaColor(object),
           pt.size = 0.5,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = T,label.size = 3, 
           title = "UMAP Plot based on ICA for 28 samples SCINA",
           do.print = T,do.return = F)
save(object, file = "data/Lung_28_20200116.Rda")
