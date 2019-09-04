########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","magrittr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
library(MAST)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)


# SLURM_ARRAY_TASK_ID
set.seed=101
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))


# samples
samples = c("combined","distal","proximal","terminal")
(con <- samples[args])

# Load Seurat
(load(file="data/Lung_24_20190824.Rda"))
Idents(object) = "orig.ident"
object %<>% subset(idents ="UNC-44-P", invert = T)
DefaultAssay(object) = "SCT"

if(con != "combined") {
    cellUse = object$conditions %in% con
    object <- object[,cellUse]
}
# FindNeighbors
npcs =100
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.6)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
TSNEPlot.1(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,
           label.repel = T,do.print = T,unique.name = "conditions",
           label.size = 4, repel = T,no.legend = T,
           title = paste("tSNE plot of",con, "samples"))

UMAPPlot.1(object, group.by="integrated_snn_res.0.6",pt.size = 1,label = T,
           label.repel = T,do.print = T,unique.name = "conditions",
           label.size = 4, repel = T,no.legend = T,
           title = paste("UMAP plot of",con, "samples"))
save(object, file = paste0("data/Lung_23",con,"_20190824.Rda"))

# Differential analysis
Idents(object) = "integrated_snn_res.0.6"
Lung_markers <- FindAllMarkers.UMI(object, logfc.threshold = 0.1, only.pos = T,
                                   test.use = "MAST")
write.csv(Lung_markers,paste0(path,"Lung_3-",con,"_markers.csv"))
Top_n = 5
top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
object %<>% ScaleData(features=unique(c(as.character(top$gene))))

DoHeatmap.1(object, marker_df = Lung_markers, Top_n = 5, do.print=T,
            angle = 0,group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
            assay = "SCT",
            label=T, cex.row=5, legend.size = 5,width=10, height=7,unique.name = "conditions",
            title = paste("Top 5 markers of each clusters in",con,"sampels"))

UMAPPlot.1(object, group.by = "conditions",split.by = "conditions",
           label = T,label.repel = T, 
           pt.size = 0.5,label.size = 4, repel = T,no.legend = T,do.print = T,
           do.return = F,title = "Cell types")