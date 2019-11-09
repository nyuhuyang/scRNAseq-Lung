########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","magrittr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
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
df_samples <- readxl::read_excel("doc/20190815_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",c(1:3,5:15)))
sample_n = intersect(sample_n, grep("COPD",df_samples$group,invert = T))
df_samples <- df_samples[sample_n,]
print(df_samples)
(samples = df_samples$sample)


#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_24_20190918.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.Seurat)

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
remove(sce_list,object_list);GC()

(remove <- which(colnames(object@meta.data) %in% "ident"))
meta.data = object@meta.data[,-remove]
object@meta.data = meta.data 
remove(meta.data);GC()
#======1.2 QC, pre-processing and normalizing the data=========================
object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = df_samples$sample)
Idents(object) = "orig.ident"
(load(file = "output/20190918/g1_24_20190918.Rda"))

object %<>% subset(subset = nFeature_RNA > 200  & #nCount_RNA > 1500 & 
                       percent.mt < 25)
# FilterCellsgenerate Vlnplot before and after filteration
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
    VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
        theme(axis.text.x = element_text(size=12),legend.position="none")
})
save(g2,file= paste0(path,"g2","_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))
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
object <- JackStraw(object, num.replicate = 20,dims = 100)
object <- ScoreJackStraw(object, dims = 1:100)

jpeg(paste0(path,"ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 85)+
    ggtitle("ElbowPlot")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18)) 
dev.off()

jpeg(paste0(path,"JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 85:95)+
    ggtitle("JackStrawPlot")+
    theme(text = element_text(size=15),	
          plot.title = element_text(hjust = 0.5,size = 18))
dev.off()
npcs =50
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8,dims.use = 1:npcs, print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs, check_duplicates = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = df_samples$sample)
Idents(object) = "orig.ident"
Idents(object) = "RNA_snn_res.0.8"
object %<>% sortIdent(numeric = T)
TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = T,
           label.repel = T,alpha = 0.9,
                 no.legend = F,label.size = 4, repel = T, title = "No Integration",
                 do.print = T)
UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = F,label.size = 4, repel = T, title = "No Integration",
                 do.print = T)
save(object, file = "data/object_orig_24_20190918.Rda")
object@assays$RNA@scale.data = matrix(0,0,0)

save(object, file = "data/Lung_24_20190918.Rda")

#======1.5 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
npcs =100
object_list %<>% lapply(function(x) {
    x %<>% RunPCA(features = object.features, verbose = FALSE)
})
options(future.globals.maxSize= object.size(object_list)*3)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features, 
                                  verbose = FALSE)
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT", 
                                  anchor.features = object.features,
                                  reference = c(1, 2), reduction = "rpca", 
                                  dims = 1:npcs)
object <- IntegrateData(anchorset = anchors, normalization.method = "SCT",dims = 1:npcs)

remove(anchors,object_list);GC()
object %<>% RunPCA(npcs = npcs, verbose = FALSE)
object <- JackStraw(object, num.replicate = 20,dims = 100)
object <- ScoreJackStraw(object, dims = 1:100)
jpeg(paste0(path,"JackStrawPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 90:100)
dev.off()
npcs =100
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
Idents(object) = "integrated_snn_res.0.8"
object %<>% sortIdent(numeric = T)
TSNEPlot.1(object, group.by="integrated_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T, alpha= 0.9,
           no.legend = F,label.size = 4, repel = T, title = "No Integration",
           do.print = T)
UMAPPlot.1(object, group.by="integrated_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T, alpha= 0.9,
                 no.legend = F,label.size = 4, repel = T, title = "No Integration",
                 do.print = T)
save(object, file = "data/Lung_24_20190918.Rda")
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

TSNEPlot.1(object = object, label = F,label.repel = T, group.by = "integrated_snn_res.0.8", 
         do.return = F, no.legend = F, title = "Integration by regions",
         pt.size = 0.3,alpha = 0.9, label.size = 5, do.print = T)

UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.0.8", 
           do.return = T, no.legend = F, title = "Integration by regions",
           pt.size = 0.2,alpha = 0.9, label.size = 5, do.print = T)

table(object$orig.ident,object$integrated_snn_res.0.6) %>% kable %>% kable_styling()
object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = "data/Lung_24_20190918.Rda")
object_data = object@assays$SCT@data
save(object_data, file = "data/Lung.data_24_20190824.Rda")


for(con in c("proximal","distal","terminal")){
    print(load(file = paste0("data/Lung_24",con,"_20190918.Rda")))
    
    meta.data = cbind.data.frame(object@meta.data,
                                 object@reductions$umap@cell.embeddings)
    meta.data = meta.data[,c("UMAP_1","UMAP_2","integrated_snn_res.0.8","manual")]
    meta.data$integrated_snn_res.0.8 = as.numeric(as.character(meta.data$integrated_snn_res.0.8))
    
    meta.data = meta.data[order(meta.data$integrated_snn_res.0.8),]
    print(colnames(meta.data))
    
    #data = as.matrix(DelayedArray::t(object@assays$SCT@data))
    #write.csv(DelayedArray::t(data), paste0(path,"object_",con,"_counts.csv"))
    write.csv(meta.data, paste0(path,"object_",con,"_meta.data.csv"))
}

load(file = paste0("data/Lung_24distal_20190918.Rda"))
DefaultAssay(object) <- 'integrated'
for(i in c(15,16)/10){
    object %<>% FindClusters(resolution = i)
    Idents(object) = paste0("integrated_snn_res.",i)
}
