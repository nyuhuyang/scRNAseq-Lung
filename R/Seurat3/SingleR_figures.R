library(Seurat)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggthemes)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file = "data/Lung_24_20190824.Rda"))
(load(file="output/singlerT_Lung_24_20190824.Rda"))

# if singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object)
if(length(singler$singler[[1]]$SingleR.single$labels) < ncol(object)){
        all.cell = colnames(object);length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = subset(object, cells = know.cell)
}

table(names(singler$singler[[1]]$SingleR.single$labels) == colnames(object))
singler$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = object@reductions$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Idents(object) # the Seurat clusters (if 'clusters' not provided)
save(singler,file="output/singlerT_Lung_24_20190824.Rda")
##############################
# add singleR label to Seurat
###############################

singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "orig.ident"  = object$orig.ident,
                       row.names = names(singler$singler[[1]]$SingleR.single$labels))
head(singlerDF)
apply(singlerDF,2,function(x) length(unique(x)))

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR::SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR::SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

table(singlerDF$singler1sub) %>% kable %>% kable_styling()
##############################
# process color scheme
##############################

apply(singlerDF[,c("singler1sub","singler1main")],2,function(x) length(unique(x)))
singlerDF[,c("singler1sub")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaData(object = object,metadata = singlerDF)
object <- AddMetaColor(object = object, label= "singler1main", colors = Singler.colors)
Idents(object) <- "singler1main"

UMAPPlot.1(object, cols = ExtractMetaColor(object),label = T, label.repel = T,pt.size = 1,
         label.size = 4, repel = T,no.legend = T,do.print = T,
         title = "Major cell types")
save(object,file="data/Lung_24_20190824.Rda")
##############################
# draw tsne plot
##############################
jpeg(paste0(path,"PlotTsne_split_lable.jpeg"), units="in", width=10, height=7,res=600)
TSNEPlot(object, cols = ExtractMetaColor(object),label = T,pt.size = 1,
         split.by = "orig.ident", label.size = 4, repel = T)+NoLegend()
dev.off()

Idents(object) <-"RNA_snn_res.0.6"
jpeg(paste0(path,"PlotTsne_split_seurat_clusters.jpeg"), units="in", width=10, height=7,res=600)
TSNEPlot(object, cols = ExtractMetaColor(object),label = T,pt.size = 1,
         group.by = "seurat_clusters",
         split.by = "orig.ident", label.size = 4, repel = T)+NoLegend()
dev.off()
