########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(kableExtra)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
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
df_samples <- readxl::read_excel("doc/20190509_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",4)))
df_samples = df_samples[sample_n,]
df_samples
(samples = df_samples$sample)


#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_8_20190807.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.Seurat)

for(i in 1:length(samples)){
    object_list[[i]]$orig.ident <- df_samples$sample[i]
    object_list[[i]]$conditions <- df_samples$conditions[i]
    object_list[[i]]$project <- df_samples$project[i]
    object_list[[i]]$tests <- df_samples$tests[i]
    Idents(object_list[[i]]) <- df_samples$sample[i]
}

#========1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_list)
object@assays$RNA@data = object@assays$RNA@data *log(2)
remove(sce_list,object_list);GC()

(remove <- which(colnames(object@meta.data) %in% "ident"))
meta.data = object@meta.data[,-remove]
object@meta.data = meta.data 

#======1.2 QC, pre-processing and normalizing the data=========================
object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = df_samples$sample)
(load(file = "output/20190807/g1_8_20190807.Rda"))

object %<>% subset(subset = nFeature_RNA > 1000  & nCount_RNA > 2000 & percent.mt < 25)
# FilterCellsgenerate Vlnplot before and after filteration
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
    VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
        theme(axis.text.x = element_text(size=15),legend.position="none")
})

save(g2,file= paste0(path,"g2_8_20190807.Rda"))
jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nFeature_RNA before filteration")+
                    scale_y_log10(limits = c(100,10000))+
                    theme(plot.title = element_text(hjust = 0.5)),
                g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                    scale_y_log10(limits = c(100,10000))+
                    theme(plot.title = element_text(hjust = 0.5))))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nCount_RNA before filteration")+
                    scale_y_log10(limits = c(500,100000))+
                    theme(plot.title = element_text(hjust = 0.5)),
                g2[[2]]+ggtitle("nCount_RNA after filteration")+ 
                    scale_y_log10(limits = c(500,100000))+
                    theme(plot.title = element_text(hjust = 0.5))))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % before filteration")+
                    ylim(c(0,50))+
                    theme(plot.title = element_text(hjust = 0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,50))+
                    theme(plot.title = element_text(hjust = 0.5))))
dev.off()

######################################

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
object@assays$RNA@scale.data = matrix(0,0,0)
object <- JackStraw(object, num.replicate = 20,dims = 100)
object <- ScoreJackStraw(object, dims = 1:85)

jpeg(paste0(path,"ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 85)+
    ggtitle("ElbowPlot")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18)) 
dev.off()

jpeg(paste0(path,"JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 60:70)+
    ggtitle("JackStrawPlot")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18)) 
dev.off()
npcs =65
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                       dims.use = 1:npcs,print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)

object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = df_samples$sample)
p0 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 no.legend = F,label.size = 4, repel = T, title = "Original")
save(p0,file= paste0(path,"p0_Original_20190807.Rda"))
#======1.5 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
options(future.globals.maxSize= 118388608000)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features, 
                                  verbose = FALSE)
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT", 
                                  anchor.features = object.features)
object <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

remove(anchors,object_list);GC()
object %<>% RunPCA(npcs = 100, verbose = FALSE)
#object <- JackStraw(object, num.replicate = 20,dims = 100)
#object <- ScoreJackStraw(object, dims = 1:100)
#jpeg(paste0(path,"JackStrawPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
#JackStrawPlot(object, dims = 90:100)
#dev.off()
npcs =65
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 1.2,
                         dims.use = 1:npcs,print.output = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
p1 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "Intergrated UMAP plot")
p2 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "Intergrated tSNE plot")
#=======1.9 summary =======================================

jpeg(paste0(path,"S1_remove_batch_tsne.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Clustering without integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p2+ggtitle("Clustering with integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

TSNEPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.0.6", 
         do.return = F, no.legend = F, title = "tSNE plot for all clusters",
         pt.size = 0.3,alpha = 1, label.size = 6, do.print = T)

UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.1.2", 
           do.return = F, no.legend = F, title = "UMAP plot for all clusters",
           pt.size = 0.2,alpha = 1, label.size = 6, do.print = T)

table(object$orig.ident,object$integrated_snn_res.0.6) %>% kable %>% kable_styling()
object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = "data/Lung_8_20190808.Rda")
object_data = object@assays$RNA@data
save(object_data, file = "data/Lung.data_8_20190808.Rda")
