########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   "magrittr","harmony"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/20200701_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
#sample_n = intersect(sample_n, grep("COPD",df_samples$group,invert = T))
print(df_samples)
(samples = df_samples$sample)
nrow(df_samples)

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_2_20200701.Rda"))
names(sce_list)
sce_list2 <- sce_list
(load(file = "data/sce_28_20200102.Rda"))
names(sce_list)
sce_list <- c(sce_list,sce_list2)
table(names(sce_list) == samples)

Seurat_list <- lapply(sce_list, as.Seurat)

for(i in 1:length(samples)){
        Seurat_list[[i]]$orig.ident <- df_samples$sample[i]
        Seurat_list[[i]]$conditions <- df_samples$conditions[i]
        Seurat_list[[i]]$group <- df_samples$group[i]
        Seurat_list[[i]]$project <- df_samples$project[i]
        Idents(Seurat_list[[i]]) <- df_samples$sample[i]
}

#========1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
object@assays$RNA@data = object@assays$RNA@data *log(2)
remove(sce_list,sce_list2,Seurat_list);GC()

(remove <- which(colnames(object@meta.data) %in% "ident"))
meta.data = object@meta.data[,-remove]
object@meta.data = object@meta.data[,-remove]
remove(meta.data);GC()
#======1.4 QC, pre-processing the data=========================
object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = df_samples$sample)
Idents(object) = "orig.ident"

object %<>% subset(subset = nFeature_RNA > 200  & #nCount_RNA > 1500 & 
                           percent.mt < 25)
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=12),legend.position="none")
})
jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                        scale_y_log10(limits = c(100,10000))+
                        theme(axis.text.x = element_text(size=10),
                              plot.title = element_text(hjust = 0.5)))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(g2[[2]]+ggtitle("nCount_RNA after filteration")+ 
                        scale_y_log10(limits = c(500,100000))+
                        theme(axis.text.x = element_text(size=10),
                              plot.title = element_text(hjust = 0.5)))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(g2[[3]]+ggtitle("mito % after filteration")+ 
                        ylim(c(0,50))+
                        theme(axis.text.x = element_text(size=10),
                              plot.title = element_text(hjust = 0.5)))
dev.off()
#======1.5 Performing SCTransform and integration =========================
Seurat_list <- SplitObject(object, split.by = "orig.ident")
# manually clean VU-29-D
#Seurat_list[["VU-29-D"]] %<>% subset(subset = nFeature_RNA > 900)

Seurat_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(Seurat_list, nfeatures = 3000)
npcs =30
Seurat_list %<>% lapply(function(x) {
        x %<>% RunPCA(features = object.features, verbose = T)
})
options(future.globals.maxSize= object.size(Seurat_list)*3)
Seurat_list <- PrepSCTIntegration(object.list = Seurat_list, anchor.features = object.features, 
                                  verbose = FALSE)

anchors <- FindIntegrationAnchors(Seurat_list, normalization.method = "SCT", 
                                  anchor.features = object.features,
                                  reference = c(1, 2), reduction = "rpca", 
                                  dims = 1:npcs)
remove(Seurat_list);GC()

object <- IntegrateData(anchorset = anchors, normalization.method = "SCT",dims = 1:npcs)
remove(anchors);GC()
# object@assays$integrated = NULL
# object %<>% subset(idents = "VU-29-D", invert = T)
# object %<>% subset(idents = c("VU-29-D","VU-35-D"), invert = T)
DefaultAssay(object) = "SCT"
object %<>% FindVariableFeatures(selection.method = "vst",
                                 num.bin = 20,
                                 mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object[["SCT"]]@misc = NULL
object[["integrated"]] = NULL
#======1.7 run UMAP =========================
object %<>% ScaleData
npcs <- 105
object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))

jpeg(paste0(path,"S1_ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = npcs)
dev.off()

object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.1,label = F,
           label.repel = F,alpha = 0.9, cols = Singler.colors,
           no.legend = T,label.size = 4, repel = T, title = "No Integration with 30 samples with 150 PCA",
           do.print = T, do.return = F, save.path = paste0(path,"Lung_30_"))
saveRDS(object, file = "data/Lung_30_20200702.rds")

object <- JackStraw(object, num.replicate = 20,dims = npcs)
object <- ScoreJackStraw(object, dims = 1:npcs)
jpeg(paste0(path,"JackStrawPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 100:110)
dev.off()

# serial resolution and generate seurat
object = readRDS(file = "data/Lung_30_20200702.rds")
DefaultAssay(object) = "SCT"
resolutions = c(seq(0.001,0.009, by = 0.001),seq(0.01,0.09, by = 0.01),seq(0.1,5, by = 0.1))
for(i in 1:length(resolutions)){
        object %<>% FindClusters(resolution = resolutions[i])
        Idents(object) = paste0("SCT_snn_res.",resolutions[i])
        UMAPPlot.1(object, group.by=paste0("SCT_snn_res.",resolutions[i]),pt.size = 0.3,label = T,
                   label.repel = T,alpha = 0.9,
                   do.return = F,
                   no.legend = T,label.size = 4, repel = T, 
                   title = paste("res =",resolutions[i]),
                   do.print = T, save.path = save.path)
        Progress(i,length(resolutions))
}
object@assays$RNA = NULL
format(object.size(object), unit = "GB")
saveRDS(object,file = "data/Lung_SCT_30_20200710.rds")

