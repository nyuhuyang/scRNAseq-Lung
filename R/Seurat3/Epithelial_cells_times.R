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
library(ggpubr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
set.seed(1001)
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

p0 <- TSNEPlot.1(Epi, group.by="RNA_snn_res.0.8",pt.size = 1,label = F,
                 no.legend = T,label.size = 4, repel = T, title = "tSNE plot before")
p1 <- UMAPPlot(Epi, group.by="RNA_snn_res.0.8",pt.size = 1,label = F,
               label.size = 4, repel = T)+ggtitle("UMAP plot before")+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))

#======1.8 RunHarmony=======================
npcs =71
jpeg(paste0(path,"RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(Epi %<>% RunHarmony.1(group.by.vars= "orig.ident", dims.use = 1:npcs,
                                   theta = NULL, plot_convergence = TRUE,#epsilon.harmony = -Inf,
                                   nclust = 100, max.iter.cluster = 100))
dev.off()

Epi %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
Epi %<>% FindClusters(reduction = "harmony",resolution = 0.3,
                       dims.use = 1:npcs,print.output = FALSE)
Epi %<>% RunTSNE(reduction = "harmony", dims = 1:npcs)
Epi %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
table(Idents(Epi))
Epi@meta.data$RNA_snn_res.0.3 <- plyr::mapvalues(Epi@meta.data$RNA_snn_res.0.3,
                                 from = 0:7,
                                 to = 1:8)
                                         
p2 <- TSNEPlot.1(Epi, group.by="RNA_snn_res.0.3",pt.size = 1,label = F,
                 label.size = 4, repel = T,title = "TSNE plot after")
p3 <- UMAPPlot(Epi, group.by="RNA_snn_res.0.3",pt.size = 1,label = F,
                 label.size = 4, repel = T)+ggtitle("UMAP plot after")+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
jpeg(paste0(path,"Epi_harmony_remove_batch.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(p0 + NoLegend(),p2 + NoLegend())
dev.off()
jpeg(paste0(path,"Epi_rerun_tsne~.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(p0 + NoLegend(),p2 + NoLegend())
dev.off()
jpeg(paste0(path,"Epi_tsne_umap~.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(p1 + NoLegend(),p3 + NoLegend())
dev.off()
#Epi$RNA_snn_res.0.8 <- plyr::mapvalues(Epi$RNA_snn_res.0.8,from = 0,to = 1)
TSNEPlot.1(Epi, group.by="RNA_snn_res.0.3",pt.size = 1,label = F,
           label.size = 4, repel = T,title = "TSNE plot res=0.3", do.print = T)
p4 <- UMAPPlot(Epi, group.by="RNA_snn_res.0.3",pt.size = 1,label = T,
               label.size = 4, repel = T)+ggtitle("UMAP plot res=0.3")+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))

jpeg(paste0(path,"Epi_umap.jpeg"), units="in", width=10, height=7,res=600)
print(p4)
dev.off()

meta.data = cbind.data.frame(Epi@meta.data,Epi@reductions$tsne@cell.embeddings,
                             Epi@reductions$umap@cell.embeddings)
colnames(meta.data)
meta.data = meta.data[,c("tSNE_1","tSNE_2","UMAP_1","UMAP_2","RNA_snn_res.0.3")]
meta.data$RNA_snn_res.0.3 = as.numeric(as.character(meta.data$RNA_snn_res.0.3))

meta.data = meta.data[order(meta.data$RNA_snn_res.0.3),]
write.csv(meta.data[,c("tSNE_1","tSNE_2","RNA_snn_res.0.3")], 
          paste0(path,"Epi_tSNE_coordinates.csv"))
write.csv(meta.data[,c("UMAP_1","UMAP_2","RNA_snn_res.0.3")], 
          paste0(path,"Epi_UMAP_coordinates.csv"))

#4 - List of markers for each cluster =================
Idents(Epi) <- "RNA_snn_res.0.3"
Epi_markers <- FindAllMarkers.UMI(Epi, logfc.threshold = 0.05,return.thresh = 0.05,
                                   only.pos = T)
table(Epi_markers$cluster)
Epi_markers = Epi_markers[Epi_markers$p_val_adj < 0.05,]
write.csv(Epi_markers,paste0(path,"Epi_markers.csv"))

Epi %<>% ScaleData(features=unique(Epi_markers$gene))
DoHeatmap.1(Epi, marker_df = Epi_markers, Top_n = 5, do.print=T, angle = 0,
            group.bar = T, title.size = 13, no.legend = F,size=5.5,hjust = 0.5,
            label= T, cex.row=5, legend.size = 5,width=10, height=7,
            title = "Top 5 markers in each clusters")
#5 - Expression data (with cluster label and sample ID for each cell; order by cluster number)
data = as.matrix(t(Epi@assays$RNA@data))
exp_data = cbind(meta.data[,4:5], data[match(rownames(meta.data),
                                              rownames(data)), ])

write.csv(t(exp_data[,-1]), paste0(path,"Epithelial_exp.csv"))
Epi@assays$RNA@scale.data = matrix(0,0,0)
save(Epi, file = "data/Epi_harmony_12_20190627.Rda")

# 1. Clusters 0 and 1 can be merged as they are barely different based on markers and both represent type 2 cells.
# 2. The following clusters can be removed: 8 – macrophages; 9 – endothelial cells; 11 – NK cells.
(load("data/Epi_harmony_12_20190624.Rda"))
table(Idents(Epi))
Epi <- subset(Epi, idents = c(8,9,11), invert = TRUE)

# 3. Once you re-do PCA, RunHarmony tSNE and UMAP, please make the list of positive markers only, with p adjusted value <0.05.

#####################
#Volcano plot - DEG analysis
####################
# FindPairMarkers
(load("data/Epi_harmony_12_20190627.Rda"))

Idents(Epi) = "RNA_snn_res.0.3"
table(Idents(Epi))
jpeg(paste0(path,"Epi_umap.jpeg"), units="in", width=10, height=7,res=600)
UMAPPlot(Epi,cols = ExtractMetaColor(Epi),label =T)+ggtitle("UMAP plot res = 0.3")
dev.off()

ident.1 <- list(c(2,3,4,6,7),
                5,
                7,
                c(3,6),
                6,
                2,
                6)
ident.2 <- list(c(1,5,8),
                1,
                c(2,3,4,6),
                c(2,4,7),
                3,
                4,
                1)
length(ident.1);length(ident.2)
subfolder <- paste0(path,"DEG/")
gde.pair <- FindPairMarkers(Epi, ident.1 = ident.1, ident.2 = ident.2,
                            logfc.threshold = 0.1, min.cells.group =3,
                            return.thresh = 0.05, only.pos = FALSE, save.path = subfolder)
gde.pair = gde.pair[gde.pair$p_val_adj< 0.05,]
write.csv(gde.pair, paste0(subfolder,"pairwise_comparision.csv"))
gde.pair = read.csv("output/20190705/Volcano_plots/DEG/pairwise_comparision.csv",row.names = 1)
head(gde.pair,10) %>% kable %>% kable_styling

titles <- c("Airway (2+3+4+6+7) vs alveolar (1+5+8)",
            "Alveolar type I cells (5) vs type II cells (1)",
            "Airway basal stem cells (7) vs airway epithelial differentiated cells (2+3+4+6)",
            "Airway secretory cells (3+6) vs other airway epithelial cells (2+4+7)",
            "Distal airway secretory cells (6) vs other airway secretory cells (3)",
            "M-ciliated cells (2) vs D-ciliated cells (4)",
            "Distal airway secretory cells (6) vs Alveolar type II cells (1)")
# Volcano plot=========
(clusters <- unique(gde.pair$cluster1.vs.cluster2))
for(i in 1:length(clusters)){
        df <- gde.pair[gde.pair$cluster1.vs.cluster2 %in% clusters[i],]
        df$log10_p_val_adj = -log10(df$p_val_adj)
        df$log10_p_val_adj[df$log10_p_val_adj == "Inf"] = 400
        left = df[df$avg_logFC < -1,]
        right = df[df$avg_logFC > 1,]
        left = rownames(left)[left$log10_p_val_adj >= head(sort(left$log10_p_val_adj,decreasing = T),15) %>% tail(1)]
        right = rownames(right)[right$log10_p_val_adj >= head(sort(right$log10_p_val_adj,decreasing = T),15) %>% tail(1)]
        g <- ggplot(df,aes(avg_logFC,log10_p_val_adj)) + 
                geom_point() + 
                ggtitle(titles[i]) + 
                ylab("-log10(p_value_adj)")+
                theme_minimal()+
                theme(plot.title = element_text(size=20,hjust = 0.5))+
                ggrepel::geom_text_repel(aes(label = gene), 
                        data=df[c(left,right),]) +
                geom_point(color = ifelse((df$avg_logFC > 1  & df$p_val_adj < 0.05) , "red",
                                          ifelse((df$avg_logFC < -1 & df$p_val_adj < 0.05), "blue","gray")))
        jpeg(paste0(path,"Volcano_plot",clusters[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
}
