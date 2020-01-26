########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","magrittr"), function(x) {
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

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/20190815_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",c(1:3,5:15)))
#sample_n = intersect(sample_n, grep("COPD",df_samples$group,invert = T))
df_samples <- df_samples[sample_n,]
print(df_samples)
(samples = df_samples$sample)
nrow(df_samples)

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_28_20200102.Rda"))
names(sce_list)
Seurat_list <- lapply(sce_list, as.Seurat)

for(i in 1:length(samples)){
    Seurat_list[[i]]$orig.ident <- df_samples$sample[i]
    Seurat_list[[i]]$conditions <- df_samples$conditions[i]
    Seurat_list[[i]]$group <- df_samples$group[i]
    Seurat_list[[i]]$project <- df_samples$project[i]
    Seurat_list[[i]]$tests <- df_samples$tests[i]
    Idents(Seurat_list[[i]]) <- df_samples$sample[i]
}

#========1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
object@assays$RNA@data = object@assays$RNA@data *log(2)
remove(sce_list,Seurat_list);GC()

(remove <- which(colnames(object@meta.data) %in% "ident"))
meta.data = object@meta.data[,-remove]
object@meta.data = object@meta.data[,-remove]
remove(meta.data);GC()
#======1.4 QC, pre-processing the data=========================
object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = df_samples$sample)
Idents(object) = "orig.ident"
(load(file = "output/20200102/g1_28_20200102.Rda"))

object %<>% subset(subset = nFeature_RNA > 200  & #nCount_RNA > 1500 & 
                       percent.mt < 25)
# FilterCellsgenerate Vlnplot before and after filteration
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
    VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
        theme(axis.text.x = element_text(size=12),legend.position="none")
})
jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nFeature_RNA before filteration")+
                    scale_y_log10(limits = c(100,10000))+
                    theme(axis.text.x = element_text(size=10),
                          plot.title = element_text(hjust = 0.5)),
                g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                    scale_y_log10(limits = c(100,10000))+
                    theme(axis.text.x = element_text(size=10),
                          plot.title = element_text(hjust = 0.5))))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nCount_RNA before filteration")+
                    scale_y_log10(limits = c(500,100000))+
                    theme(axis.text.x = element_text(size=10),
                          plot.title = element_text(hjust = 0.5)),
                g2[[2]]+ggtitle("nCount_RNA after filteration")+ 
                    scale_y_log10(limits = c(500,100000))+
                    theme(axis.text.x = element_text(size=10),
                          plot.title = element_text(hjust = 0.5))))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % before filteration")+
                    ylim(c(0,50))+
                    theme(axis.text.x = element_text(size=10),
                          plot.title = element_text(hjust = 0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,50))+
                    theme(axis.text.x = element_text(size=10),
                          plot.title = element_text(hjust = 0.5))))
dev.off()
save(g2,file= paste0(path,"g2","_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))


# After removing unwanted cells from the dataset, the next step is to normalize the data.
#object <- NormalizeData(object = object, normalization.method = "LogNormalize", 
#                      scale.factor = 10000)
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20,nfeatures = 3000,
                               mean.cutoff = c(0.1, 8), 
                               dispersion.cutoff = c(1, Inf))

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

set.seed(100)
Seurat_list <- SplitObject(object, split.by = "orig.ident")

QC.list = read.csv(paste0(path,"test28_QC_list.csv"),row.names = 1, stringsAsFactors = F)
args = ""
if(is.na(args[2])){
    message("QC")
    nCount_RNA <- sapply(Seurat_list, function(x) mean(x$nCount_RNA))
    nFeature_RNA <- sapply(Seurat_list, function(x) mean(x$nFeature_RNA))
    
    QC.list <- cbind(QC.list, nCount_RNA, nFeature_RNA,
                     row.names = df_samples$sample)
    write.csv(QC.list,paste0(path,"test28_QC_list.csv"))
    QC.list %>% kable() %>% kable_styling()
    remove(QC.list,nCount_RNA,nFeature_RNA);GC()
}

remove(object);GC()

#======1.5 Performing SCTransform and integration =========================
Seurat_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(Seurat_list, nfeatures = 3000)
npcs =30
Seurat_list %<>% lapply(function(x) {
    x %<>% RunPCA(features = object.features, verbose = FALSE)
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

npcs = 100
object %<>% RunPCA(npcs = npcs, verbose = FALSE)
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
Idents(object) = "integrated_snn_res.0.8"
object %<>% sortIdent(numeric = T)
p2 <- TSNEPlot.1(object, group.by="integrated_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T, alpha= 0.9,
           no.legend = F,label.size = 4, repel = T, title = "With Integration",
           do.print = T, do.return = T)
p3 <- UMAPPlot.1(object, group.by="integrated_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T, alpha= 0.9,
           no.legend = F,label.size = 4, repel = T, title = "With Integration",
           do.print = T, do.return = T)
g3 <- UMAPPlot.1(object, group.by="conditions",pt.size = 1,label = T,
                 cols = c('#601b3f','#3D9970','#FF4136','#FF851B'),
                 label.repel = T, alpha= 0.9,
                 no.legend = F,label.size = 4, repel = T, title = "With Integration",
                 do.print = T, do.return = T)
save(object, file = "data/Lung_28_20200102.Rda")


#======1.6 without intergration =========================
DefaultAssay(object) = "SCT"
#object <- FindVariableFeatures(object = object, selection.method = "vst",
#                               num.bin = 20,
#                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
VariableFeatures(object) = rownames(object@assays$SCT@scale.data)
#object <- ScaleData(object = object,features = rownames(object))
object %<>% RunICA(verbose =F,nics = 100)
object <- JackStraw(object, reduction = "ica",num.replicate = 20,dims = 100)
object <- ScoreJackStraw(object, dims = 1:100)
jpeg(paste0(path,"JackStrawPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 90:100)
dev.off()
npcs =100
npcs =100
object %<>% FindNeighbors(reduction = "ica",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunTSNE(reduction = "ica", dims = 1:npcs, check_duplicates = FALSE)
object %<>% RunUMAP(reduction = "ica", dims = 1:npcs)

object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = df_samples$sample)
Idents(object) = "RNA_snn_res.0.8"
object %<>% sortIdent(numeric = T)
p0 <- TSNEPlot.1(object, group.by="integrated_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T,alpha = 0.9,
             no.legend = F,label.size = 4, repel = T, title = "No Integration",
             do.print = T, do.return = T)
p1 <- UMAPPlot.1(object, group.by="integrated_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = F,label.size = 4, repel = T, title = "No Integration",
           do.print = T, do.return = T)

g1 <- UMAPPlot.1(object, group.by="conditions",pt.size = 1,label = T,
                 cols = c('#601b3f','#3D9970','#FF4136','#FF851B'),
                 label.repel = T,alpha = 0.9,
                 no.legend = F,label.size = 4, repel = T, title = "No Integration",
                 do.print = T, do.return = T)
save(object, file = "data/Lung_28_20200103.Rda")

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

jpeg(paste0(path,"S1_conditions_umap.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1+ggtitle("Clustering without integration")+NoLegend()+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          g3+ggtitle("Clustering with integration")+NoLegend()+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()





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

(load(file = "data/Lung_28_20200116.Rda"))
DefaultAssay(object) <- 'SCT'
object@reductions= readRDS(file = "output/20200114/reductions_4_SCT_2000.rds")
res = c(6:20)/10
for(i in seq_along(res)){
    object %<>% FindClusters(resolution = res[i])
    Idents(object) = paste0("RNA_snn_res.",res[i])
    TSNEPlot.1(object, group.by=paste0("SCT_snn_res.",res[i]),pt.size = 0.3,label = T,
               label.repel = T,alpha = 0.9,cols = c(Singler.colors,Singler.colors),
               do.return = F,
               no.legend = T,label.size = 4, repel = T, 
               title = paste("res =",res[i]," based on ICA"),
               do.print = T)
    Progress(i,length(res))
}
