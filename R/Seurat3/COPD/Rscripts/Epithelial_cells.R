########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(gplots)
library(MAST)
library(ggpubr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
set.seed(1001)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# samples
samples = c("proximal","distal","terminal")
(con <- samples[args])

# 5.0 Preliminaries: Load the data
(load(file = paste0("data/Lung_23",con,"_20190824.Rda")))
Idents(object) = "cell.type"
object %<>% subset(idents=c("CD4+ T-cells", "CD8+ T-cells", "DC"), invert = TRUE)
Idents(object) = "integrated_snn_res.0.6"
Epi <- subset(object, idents=c(2, 7, 8, 10, 13, 20, 9))
remove(object);GC()
cluster_9 <- subset(Epi, idents = 9)
DefaultAssay(cluster_9) = "RNA"
cluster_9 <- subset(cluster_9, subset = (KRT19 == 0 | KRT5 == 0))
Epi.cells <- colnames(Epi)
Epi.cells = Epi.cells[!(Epi.cells %in% colnames(cluster_9))]
Epi <- subset(Epi, cells = Epi.cells)
#======1.5 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(Epi, split.by = "orig.ident")
remove(Epi);GC()
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
options(future.globals.maxSize= object.size(object_list)*1.5)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features, 
                                  verbose = FALSE)
k.filter <- min(200, min(sapply(object_list, ncol)))-1
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT", 
                                  anchor.features = object.features, k.filter = k.filter)
object <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

remove(anchors,object_list);GC()
object %<>% RunPCA(npcs = 100, verbose = FALSE)
npcs =100
object %<>% FindNeighbors
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
for(i in 4:14/10) {
        object %<>% FindClusters(resolution = i)
        UMAPPlot.1(object = object, label = T,label.repel = T, group.by = paste0("integrated_snn_res.",i), 
                   do.return = F, no.legend = T, title = paste("UMAP plot of",con,"samples"),
                   pt.size = 0.2,alpha = 1, label.size = 5, do.print = T, unique.name = "conditions")
}
TSNEPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.0.8", 
           do.return = F, no.legend = T, title = paste("tSNE plot of",con,"samples"),
           pt.size = 0.3,alpha = 1, label.size = 5, do.print = T, unique.name = "conditions")
Epi <- object
save(Epi, file = paste0("data/Epi_23",con,"_20190904.Rda"))

Idents(Epi) = "integrated_snn_res.0.8"
Lung_markers <- FindAllMarkers.UMI(Epi, logfc.threshold = 0.1, only.pos = T,
                                   test.use = "MAST")
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
write.csv(Lung_markers,paste0(path,"Epi_23-",con,"_markers.csv"))

Top_n = 5
top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
Epi %<>% ScaleData(features=unique(c(as.character(top$gene))))

DoHeatmap.1(Epi, marker_df = Lung_markers, Top_n = 5, do.print=T,
            angle = 0,group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
            assay = "SCT",
            label=T, cex.row=5, legend.size = 5,width=10, height=7,unique.name = "conditions",
            title = paste("Top 5 markers of each clusters in",con,"sampels"))