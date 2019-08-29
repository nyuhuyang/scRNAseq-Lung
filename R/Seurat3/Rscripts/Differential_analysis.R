########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
library(MAST)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
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

# 3.1.1 load data
(load(file="data/Lung_24_20190824.Rda"))
Idents(object) = "group"
object %<>% subset(idents ="UNC-44", invert = T)
DefaultAssay(object) = "SCT"
Idents(object) = "integrated_snn_res.0.6"

if(con == "combined") {
    sub_object <- object
} else {
    cellUse = object$conditions %in% con
    sub_object <- object[,cellUse]
}

p1 <- UMAPPlot.1(object = sub_object, label = F,label.repel = F, group.by = "integrated_snn_res.0.6",cols = ExtractMetaColor(object),do.return = T, no.legend = F, title = paste("UMAP plot for all clusters in",con),pt.size = 0.2,alpha = 1, label.size = 5, do.print = F,unique.name = T)
jpeg(paste0(path,"UMAPPlot_clusters_",con,"~.jpeg"),units='in', width=10, height=7,res=600)
print(p1)
dev.off()

p2 <- TSNEPlot.1(object = sub_object, label = F,label.repel = F, group.by = "integrated_snn_res.0.6",cols = ExtractMetaColor(object),do.return = T, no.legend = F, title = paste("tSNE plot for all clusters in",con),pt.size = 0.2,alpha = 1, label.size = 5, do.print = F,unique.name = T)
jpeg(paste0(path,"TSNEPlot_clusters_",con,"~.jpeg"),units='in', width=10, height=7,res=600)
print(p2)
dev.off()

Idents(sub_object) = "integrated_snn_res.0.6"
Lung_markers <- FindAllMarkers.UMI(sub_object, logfc.threshold = 0.75, only.pos = T,
                                   test.use = "MAST")
write.csv(Lung_markers,paste0(path,"Lung_24-",con,"_markers.csv"))

Top_n = 5
top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
sub_object %<>% ScaleData(features=unique(c(as.character(top$gene))))

DoHeatmap.1(sub_object, marker_df = Lung_markers, Top_n = 5, do.print=F,
            angle = 0,group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
            assay = "SCT",
            label=T, cex.row=5, legend.size = 5,width=10, height=7,unique.name = "conditions",
            title = paste("Top 5 markers of each clusters in",con,"sampels"))
