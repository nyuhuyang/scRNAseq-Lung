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

(load(file = paste0("data/Lung_16_distal_20191022.Rda")))
print(unique(object@meta.data$conditions))

unique(gsub("_.*","",colnames(object)))
# - Table: number of cells per cluster (per each sample and total)
df <- table(object$labels, object$orig.ident) %>% 
        as.data.frame()
colnames(df) = c("labels","samples","Freq")
df %<>% spread("samples","Freq")
rownames(df) = df$labels
df = df[,-1]
write.csv(df, paste0(path,"Lung_16_distal_cell.types_by_samples.csv"))

#  UMAP coordinates &  Expression data =============
Idents(object) = "conditions"
distal  = subset(object, idents = "distal")
COPD  = subset(object, idents = "COPD")
for(x in list(distal,COPD)){
        meta.data = cbind.data.frame(x@meta.data,
                                     x@reductions$umap@cell.embeddings)
        meta.data = meta.data[,c("UMAP_1","UMAP_2","RNA_snn_res.0.8","labels")]
        meta.data$RNA_snn_res.0.8 = as.numeric(as.character(meta.data$RNA_snn_res.0.8))
        
        meta.data = meta.data[order(meta.data$RNA_snn_res.0.8),]
        print(colnames(meta.data))
        
        data = as.matrix(DelayedArray::t(x@assays$RNA@data))
        write.csv(DelayedArray::t(data), paste0(path,"Lung_16_",unique(x$conditions),"_counts.csv"))
        write.csv(meta.data, paste0(path,"Lung_16_",unique(x$conditions),"_meta.data.csv"))
}


# - Table: number of expressed genes per cluster (per each sample and total) ==============
df <- table(object$labels, object$orig.ident) %>% as.data.frame()
colnames(df) =c("labels","samples","nGene")

df$labels <- as.character(df$labels)
df$samples <- as.character(df$samples)

for(i in seq_len(nrow(df))) {
        cells <- object$labels %in% df[i,"labels"] & object$orig.ident %in% df[i,"samples"]
        df[i,"nGene"] = as.integer(mean(object@meta.data[cells,"nFeature_RNA"]))
        svMisc::progress(i/nrow(df)*100)
}
df %<>% spread("samples","nGene")
df[is.na(df)] = 0
write.csv(df,paste0(path,"Lung_16_distal_nGene_by_samples.csv"))

