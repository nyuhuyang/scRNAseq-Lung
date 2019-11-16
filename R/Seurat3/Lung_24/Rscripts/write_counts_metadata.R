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
library(Matrix)
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

samples = c("proximal","distal","terminal")
(con <- samples[args])

path <- paste0(path,con,"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#(load(file = paste0("data/Lung_24-",con,"_20191004.Rda")))
(load(file = paste0("data/Lung_24_20190918.Rda")))
print(unique(object@meta.data$conditions))

# - Table: number of cells per cell types (per each sample and total)
df <- table(object$cell.types, object$orig.ident) %>% 
        as.data.frame()
colnames(df) = c("cell.types","samples","Freq")
df %<>% spread("samples","Freq")
rownames(df) = df$cell.types
df = df[,-1]
#write.csv(df, paste0(path,"Lung_24-",con,"_cell.types_by_samples.csv"))
write.csv(df, paste0(path,"Lung_24_cell.types_by_samples.csv"))

# - Table: number of cells per cluster (per each sample and total)
object$RNA_snn_res.0.8 %<>% as.numeric()
df <- table(object$RNA_snn_res.0.8, object$orig.ident) %>% 
        as.data.frame()
colnames(df) = c("RNA_snn_res.0.8","samples","Freq")
df %<>% spread("samples","Freq")
rownames(df) = df$RNA_snn_res.0.8
df = df[,-1]
#write.csv(df, paste0(path,"Lung_24-",con,"_cell.types_by_samples.csv"))
write.csv(df, paste0(path,"Lung_24_clusters_by_samples.csv"))

#  UMAP, tSNE coordinates 
meta.data = cbind.data.frame(object@meta.data,
                             object@reductions$umap@cell.embeddings,
                             object@reductions$tsne@cell.embeddings)
meta.data = meta.data[,c("UMAP_1","UMAP_2","tSNE_1","tSNE_2","RNA_snn_res.0.8","cell.types")]
meta.data$RNA_snn_res.0.8 = as.numeric(as.character(meta.data$RNA_snn_res.0.8))

meta.data = meta.data[order(meta.data$RNA_snn_res.0.8),]
print(colnames(meta.data))
write.csv(meta.data[,c("UMAP_1","UMAP_2","RNA_snn_res.0.8","cell.types")], 
          paste0(path,"Lung_24_UMAP_coordinates.csv"))
write.csv(meta.data[,c("tSNE_1","tSNE_2","RNA_snn_res.0.8","cell.types")], 
          paste0(path,"Lung_24_tSNE_coordinates.csv"))
# Expression data =============

data = summary(object@assays$RNA@data)
system.time(df <- spread(data, key = j,value = x,drop=FALSE,
                         fill = 0L))
fwrite(df, file= paste0(path,"Lung_counts.csv"))

# - Table: number of expressed genes per cluster (per each sample and total) ==============
df <- table(object$RNA_snn_res.0.8, object$orig.ident) %>% as.data.frame()
colnames(df) =c("cluster_res.0.8","samples","nGene")

df$cluster_res.0.8 <- as.character(df$cluster_res.0.8)
df$samples <- as.character(df$samples)

for(i in seq_len(nrow(df))) {
        cells <- object$RNA_snn_res.0.8 %in% df[i,"cluster_res.0.8"] & object$orig.ident %in% df[i,"samples"]
        df[i,"nGene"] = as.integer(mean(object@meta.data[cells,"nFeature_RNA"]))
        svMisc::progress(i/nrow(df)*100)
}
df %<>% spread("samples","nGene")
df[is.na(df)] = 0
#write.csv(df,paste0(path,"Lung_24-",con,"_nGene_by_samples.csv"))
write.csv(df,paste0(path,"Lung_24_nGene_by_samples.csv"))

