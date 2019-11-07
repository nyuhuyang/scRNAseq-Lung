########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(magrittr)
library(DoubletFinder)
library(kableExtra)
library(gplots)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# clustering dendrogram (Spearman or Pearson) for all cell types/clusters (16 D+COPD) 
# based on top (50-100?) genes in each clust
(load(file = "data/Lung_16_distal_20191022.Rda"))
Lung_markers = read.csv(file = "Yang/distal_COPD/DE analysis/Lung_16_RNA_snn_res.0.8_markers.csv",
                        row.names = 1,stringsAsFactors = F)
Lung_markers = read.csv(file = "Yang/distal_COPD/DE analysis/Lung_16_distal_labels_markers.csv",
                        row.names = 1,stringsAsFactors = F)
Lung_markers1 = read.csv(file = "Yang/distal_COPD/DE analysis/Lung_16_COPD_labels_markers.csv",
                        row.names = 1,stringsAsFactors = F)
Lung_markers = rbind(Lung_markers,Lung_markers1)
#Idents(object) = "RNA_snn_res.0.8"
Idents(object) = "labels"
object %<>% sortIdent()
table(Idents(object))
object_exp <- AverageExpression(object,assays = "RNA")
Top_n = 500
top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
exp =  object_exp$RNA[unique(c(as.character(top$gene))),]
exp = exp %>% t %>% scale %>% t

features = c(as.character(top$gene))
featuresNum <- make.unique(features, sep = ".")
exp %<>% MakeUniqueGenes(features = features)

hc <- hclust(as.dist(1-cor(exp, method="spearman")), method="complete")
cc = as.character(as.numeric(as.factor(hc$labels)))

jpeg(paste0(path,"Heatmap2_celltype_",Top_n,".jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(exp),
          Colv = as.dendrogram(hc),
          Rowv= TRUE,
          ColSideColors = cc, 
          trace ="none",
          dendrogram = "column",
          key.xlab = "scale log nUMI",
          cexRow = 0.5,
          margins = c(12,5),
          breaks = seq(-3,3,length.out = 101),
          col = bluered,
          main = paste("Clustering dendrogram for all cell types\n (16 D+COPD) based on top",Top_n, 
                       "genes"))
dev.off()
