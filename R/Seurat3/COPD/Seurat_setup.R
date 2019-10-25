########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","magrittr",
                   "BiocSingular"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/20190815_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",c(1:3,5:15)))
sample_n = intersect(sample_n, grep("distal|COPD",df_samples$conditions))
df_samples <- df_samples[sample_n,]
print(df_samples)
(attach(df_samples))
(samples <- df_samples$sample)


#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_4_20191009.Rda"))
names(sce_list)
sce_list_COPD <- sce_list
(load(file = "data/sce_24_20191016.Rda"))
names(sce_list)
sce_list_all <- c(sce_list_COPD,sce_list)
sce_list <- sce_list_all[samples]
names(sce_list)

object_list <- lapply(sce_list, as.Seurat)
remove(sce_list, sce_list_all,sce_list_COPD);GC()
for(i in 1:length(samples)){
    object_list[[i]]$orig.ident <- df_samples$sample[i]
    object_list[[i]]$conditions <- df_samples$conditions[i]
    object_list[[i]]$group <- df_samples$group[i]
    object_list[[i]]$project <- df_samples$project[i]
    object_list[[i]]$tests <- df_samples$tests[i]
    Idents(object_list[[i]]) <- df_samples$sample[i]
}
#========1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_list)
object@assays$RNA@data = object@assays$RNA@data *log(2)
remove(object_list);GC()

(remove <- which(colnames(object@meta.data) %in% "ident"))
meta.data = object@meta.data[,-remove]
object@meta.data = meta.data 
remove(meta.data);GC()
#======1.2 QC, pre-processing and normalizing the data=========================
Idents(object) = "orig.ident"

# After removing unwanted cells from the dataset, the next step is to normalize the data.
#object <- NormalizeData(object = object, normalization.method = "LogNormalize", 
#                      scale.factor = 10000)
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(object), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
jpeg(paste0(path,"VariableFeaturePlot.jpeg"), units="in", width=10, height=7,res=600)
print(plot2)
dev.off()
hvf.info <- HVFInfo(object = object)
hvf.info = hvf.info[VariableFeatures(object),]
write.csv(hvf.info, file = paste0(path,"high_variable_genes.csv"))
#======1.3 1st run of pca-tsne  =========================
object <- ScaleData(object = object,features = rownames(object))
object <- RunPCA(object, features = VariableFeatures(object),verbose =F,npcs = 100)

npcs =50
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8,dims.use = 1:npcs, print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs, check_duplicates = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = df_samples$sample)
p0 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 no.legend = F,label.size = 4, repel = T, title = "Original")
p1 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 no.legend = F,label.size = 4, repel = T, title = "Original")
save(object, file = "data/Lung_16_distal_20191022.Rda")

#======1.5 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
#remove(object);GC()
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
options(future.globals.maxSize= object.size(object_list)*1.5)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features, 
                                  verbose = FALSE)
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT", 
                                  anchor.features = object.features)
object <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

remove(anchors,object_list);GC()
object %<>% RunPCA(npcs = 100, verbose = FALSE)
#object <- JackStraw(object, num.replicate = 20,dims = 100)
#object <- ScoreJackStraw(object, dims = 1:100)

#jpeg(paste0(path,"JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
#JackStrawPlot(object, dims = 80:90)+
#    ggtitle("JackStrawPlot")+
#    theme(text = element_text(size=15),	
#          plot.title = element_text(hjust = 0.5,size = 18))
#dev.off()
library(RcppParallel)
setThreadOptions(numThreads=1)
npcs = 80
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
p2 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "Intergrated tSNE plot")
p3 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "Intergrated UMAP plot")
object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = "data/Lung_16_distal_20191017.Rda")
#=======1.9 summary =======================================
jpeg(paste0(path,"S1_remove_batch_tsne.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Clustering without integration")+NoLegend()+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p2+ggtitle("Clustering with integration")+NoLegend()+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

jpeg(paste0(path,"S1_remove_batch_umap.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1+ggtitle("Clustering without integration")+NoLegend()+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p3+ggtitle("Clustering with integration")+NoLegend()+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

TSNEPlot.1(object = object, label = T,label.repel = T, 
         do.return = F, no.legend = F, title = "tSNE plot for all clusters",
         pt.size = 0.3,alpha = 1, label.size = 5, do.print = T)

UMAPPlot.1(object = object, label = T,label.repel = T,cols = Singler.colors,
           do.return = T, no.legend = F, title = "UMAP plot for all clusters",
           pt.size = 0.2,alpha = 1, label.size = 5, do.print = T)

object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = "data/Lung_16_distal_20191017.Rda")

DefaultAssay(object) <- 'RNA'
for(i in c(4:16)/10){
    object %<>% FindClusters(resolution = i)
    Idents(object) = paste0("RNA_snn_res.",i)
    UMAPPlot.1(object = object, label = T,label.repel = T, 
               group.by = paste0("RNA_snn_res.",i), 
               cols = Singler.colors,legend.size = 15,
               do.return = F, no.legend = T, title = paste("UMAP plot for all clusters in resolution",i),
               pt.size = 0.2,alpha = 0.85, label.size = 3, do.print = T,
               unique.name = "conditions")
}

