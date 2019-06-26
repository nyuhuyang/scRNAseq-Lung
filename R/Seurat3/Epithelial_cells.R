#1 - Re-run tSNE separately for epithelial cells (KRT19+, EPCAM+): clusters 0, 9, 11, 13, 15, 18, 25, 26. If possible, not just remove other clusters from existing tSNE, but re-run the analysis for cells belonging to epithelial clusters, so we may possibly see different distribution of clusters on tSNE map. I suggest to exclude cluster 14 from this analysis, since it has high expression of macrophage markers.
#2 - tSNE coordinates (with cluster label and sample ID columns; ordered by cluster number)
#3 - Variance plot and gene list for this specific analysis
#4 - List of markers for each cluster
#5 - Expression data (with cluster label and sample ID for each cell; order by cluster number)
########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
library(MAST)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# 3.1.1 load data
# Rename ident
(load(file = "data/Lung_harmony_12_20190614.Rda"))
Idents(object) <- "RNA_snn_res.1.2"

Epi <- subset(object, idents = c(0,9,11,13,15,18,25,26))
remove(object);GC()

p0 <- TSNEPlot.1(Epi, group.by="RNA_snn_res.1.2",pt.size = 1,label = T,
                 label.size = 4, repel = T,title = "Before Re-run tSNE")
# 3.1.2
Epi %<>% FindVariableFeatures(selection.method = "vst",
                               num.bin = 20,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(Epi), 20)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Epi)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
jpeg(paste0(path,"Epi_VariableFeaturePlot.jpeg"), units="in", width=10, height=7,res=600)
print(plot2)
dev.off()
hvf.info <- HVFInfo(object = Epi)
hvf.info = hvf.info[VariableFeatures(Epi),]
write.csv(hvf.info, file = paste0(path,"Epi_high_variable_genes.csv"))

# Re-run tsne plot
#======1.3 1st run of pca-tsne  =========================

Epi %<>% ScaleData(features = rownames(Epi))
Epi %<>% RunPCA(features = VariableFeatures(Epi),verbose =F,npcs = 100)
Epi %<>% JackStraw(num.replicate = 20,dims = 100)
Epi %<>% ScoreJackStraw(dims = 1:85)

jpeg(paste0(path,"Epi_JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(Epi, dims = 68:78)+
        ggtitle("JackStrawPlot for Epithelial cells")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18)) 
dev.off()

npcs =71
Epi %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
Epi %<>% FindClusters(reduction = "pca",resolution = 0.8,
                       dims.use = 1:npcs,print.output = FALSE)
Epi %<>% RunTSNE(reduction = "pca", dims = 1:npcs)

p1 <- TSNEPlot.1(Epi, group.by="orig.ident",pt.size = 1,label = F,
                 no.legend = T,label.size = 4, repel = T, title = "Original")
#======1.8 RunHarmony=======================
jpeg(paste0(path,"RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(Epi %<>% RunHarmony.1(group.by.vars= "orig.ident", dims.use = 1:npcs,
                                   theta = NULL, plot_convergence = TRUE,#epsilon.harmony = -Inf,
                                   nclust = 100, max.iter.cluster = 100))
dev.off()

Epi %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
Epi %<>% FindClusters(reduction = "harmony",resolution = 0.8,
                       dims.use = 1:npcs,print.output = FALSE)
Epi %<>% RunTSNE(reduction = "harmony", dims = 1:npcs)
Epi %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)

p2 <- TSNEPlot.1(Epi, group.by="RNA_snn_res.0.8",pt.size = 1,label = T,
                 label.size = 4, repel = T,title = "TSNE plot")
p3 <- UMAPPlot(Epi, group.by="RNA_snn_res.0.8",pt.size = 1,label = T,
                 label.size = 4, repel = T)+ggtitle("UMAP plot")+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))

jpeg(paste0(path,"Epi_harmony_remove_batch.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(p1 + NoLegend(),p2 + NoLegend())
dev.off()
jpeg(paste0(path,"Epi_rerun_tsne.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(p0 + NoLegend(),p2 + NoLegend())
dev.off()
jpeg(paste0(path,"Epi_tsne_umap.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(p2 + NoLegend(),p3 + NoLegend())
dev.off()
meta.data = cbind.data.frame(Epi@meta.data,Epi@reductions$tsne@cell.embeddings,
                             Epi@reductions$umap@cell.embeddings)
colnames(meta.data)
meta.data = meta.data[,c("tSNE_1","tSNE_2","UMAP_1","UMAP_2","RNA_snn_res.0.8")]
meta.data$RNA_snn_res.0.8 = as.numeric(as.character(meta.data$RNA_snn_res.0.8))

meta.data = meta.data[order(meta.data$RNA_snn_res.0.8),]
write.csv(meta.data[,c("tSNE_1","tSNE_2","RNA_snn_res.0.8")], 
          paste0(path,"Epi_tSNE_coordinates.csv"))
write.csv(meta.data[,c("UMAP_1","UMAP_2","RNA_snn_res.0.8")], 
          paste0(path,"Epi_UMAP_coordinates.csv"))

#4 - List of markers for each cluster =================
Idents(Epi) <- "RNA_snn_res.0.8"
Epi_markers <- FindAllMarkers.UMI(Epi, logfc.threshold = 0.1,
                                   only.pos = F)
table(Epi_markers$cluster)
write.csv(Epi_markers,paste0(path,"Epi_markers.csv"))

Epi %<>% ScaleData(features=unique(Epi_markers$gene))
DoHeatmap.1(Epi, marker_df = Epi_markers, Top_n = 5, do.print=T, angle = 0,
            group.bar = T, title.size = 13, no.legend = F,size=0,hjust = 0.5,
            label=F, cex.row=5, legend.size = 5,width=10, height=7,
            title = "Top 5 markers in each clusters")
#5 - Expression data (with cluster label and sample ID for each cell; order by cluster number)
data = as.matrix(t(Epi@assays$RNA@data))
exp_data = cbind(meta.data[,4:5], data[match(rownames(meta.data),
                                              rownames(data)), ])

exp_data[1:4,1:5]
write.csv(t(exp_data[,-1]), paste0(path,"Epithelial_exp.csv"))
Epi@assays$RNA@scale.data = matrix(0,0,0)
save(Epi, file = "data/Epi_harmony_12_20190624.Rda")
