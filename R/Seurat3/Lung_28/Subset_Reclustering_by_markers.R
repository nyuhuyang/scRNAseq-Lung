########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot",
                   "magrittr","harmony"), function(x) {
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
#======1.1 read marker genes =========================

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
sheets = c("T20","IE", "DS", "3C")
s = 2
df_samples <- readxl::read_excel("doc/Cell type markers for UMAP re-clustering.xlsx",
                                 sheet = sheets[s])

save.path <- paste0(path,sheets[s],"/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

#======1.2 load  Seurat =========================
(load(file = "data/Lung_28_20200102.Rda"))
object[["integrated"]] = NULL
object@neighbors = list()
object@reductions = list()
DefaultAssay(object) = "SCT"
# After removing unwanted cells from the dataset, the next step is to normalize the data.
#object <- NormalizeData(object = object, normalization.method = "LogNormalize", 
#                      scale.factor = 10000)
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20,nfeatures = 10000,
                               mean.cutoff = c(0.1, 8), 
                               dispersion.cutoff = c(1, Inf))
# keep the VariableFeatures in marker list only
marker_genes = FilterGenes(object, df_samples$genes)
VariableFeatures(object) %<>% .[. %in% marker_genes]
length(VariableFeatures(object))
# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(object), 20)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
jpeg(paste0(save.path,"VariableFeaturePlot.jpeg"), units="in", width=10, height=7,res=600)
print(plot2)
dev.off()
hvf.info <- HVFInfo(object = object)
hvf.info = hvf.info[VariableFeatures(object),]
write.csv(hvf.info, file = paste0(save.path,"high_variable_genes.csv"))

#======1.7 harmony =========================
object %<>% ScaleData
object %<>% RunPCA(verbose = T,npcs = 100, features = VariableFeatures(object))
PCAPlot.1(object, do.print = T, save.path = save.path)

npcs= 100
object %<>% JackStraw(num.replicate = 20,dims = 100)
object %<>% ScoreJackStraw(dims = 1:npcs)
a <- seq(1,100, by = 10)
b <- a+9
for(i in seq_along(a)){
    jpeg(paste0(save.path,"JackStrawPlot_",i,"_", a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a[i]:min(b[i],100)))
    Progress(i,length(a))
    dev.off()
}
saveRDS(object, file = paste0("data/Lung_28_",args,"-",g_name,"_20200206.rds"))

npcs = 100
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
theta = 2, plot_convergence = TRUE,
nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs, check_duplicates = FALSE)
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)

lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun) 
    fun(object, group.by="orig.ident",pt.size = 0.5,label = F,
        cols = ExtractMetaColor(object),
        label.repel = T,alpha = 0.9,
        no.legend = T,label.size = 4, repel = T, title = "Harmony Integration",
        do.print = T, do.return = F,save.path = save.path))

meta.data = readRDS(file = "output/20200131/cell_types.rds")
meta.data = meta.data[rownames(object@meta.data),]
object[["cell_types"]] = meta.data$cell_types
object[["cell_types.colors"]] = meta.data$cell_types.colors

Idents(object) = "cell_types"
lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
    fun(object, group.by="cell_types",pt.size = 0.5,label = T,
        label.repel = T,alpha = 0.9,cols = ExtractMetaColor(object),
        no.legend = T,label.size = 4, repel = T, title = "Harmony Integration",
        do.print = T, do.return = F,save.path = save.path))

lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
    fun(object, group.by="conditions",pt.size = 0.5,label = T,
        cols = c('#601b3f','#3D9970','#FF4136','#FF851B'),
        label.repel = T, alpha= 0.9,
        no.legend = F,label.size = 4, repel = T, title = "Harmony Integration",
        do.print = T, do.return = F,save.path = save.path))

object@assays$RNA@scale.data = matrix(0,0,0)
object@assays$SCT@scale.data = matrix(0,0,0)
object@assays$integrated@scale.data = matrix(0,0,0)
format(object.size(object[["RNA"]]),unit = "GB")
save(object, file = "data/Lung_28_harmony_rmD_20200205.Rda")
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

TSNEPlot.1(object = object, label = T,label.repel = T, group.by = "orig.ident",
           do.return = F, no.legend = F, title = "Harmony intergrated tSNE plot",
           pt.size = 0.3,alpha = 1, label.size = 5, do.print = T)

UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "SCT_snn_res.0.8",
           do.return = F, no.legend = F, title = "Harmony intergrated UMAP plot",
           pt.size = 0.2,alpha = 1, label.size = 5, do.print = T)


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

(load(file = "data/Lung_28_harmony_20200131.Rda"))
DefaultAssay(object) <- 'SCT'
options(scipen=999)
(res = c(seq(0.0001,0.0009,by = 0.0001)))
    
for(i in seq_along(res)){
    object %<>% FindClusters(resolution = res[i])
    Idents(object) = paste0("SCT_snn_res.",res[i])
    
    lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
        fun(object, group.by=paste0("SCT_snn_res.",res[i]),pt.size = 0.3,label = T,
               label.repel = T,alpha = 0.9,cols = c(Singler.colors,Singler.colors),
               do.return = F,
               no.legend = T,label.size = 4, repel = T, 
               title = paste("res =",res[i]," based on harmony"),
               do.print = T))
    Progress(i,length(res))
}

format(object.size(object@assays$SCT), units = "GB")
