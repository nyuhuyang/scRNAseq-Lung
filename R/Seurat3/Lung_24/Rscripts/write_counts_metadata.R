########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","magrittr","MAST"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
# SLURM_ARRAY_TASK_ID
set.seed=101
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

regions = c("proximal","distal","terminal")
(con <- regions[args])

path <- paste0(path,con,"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#(load(file = paste0("data/Lung_24-",con,"_20191004.Rda")))
(load(file = paste0("data/Lung_24_20191128.Rda")))
print(unique(object@meta.data$conditions))

# - Table: number of cells per cell types (per each sample and total)
df <- table(object$cell.types, object$orig.ident) %>% 
        as.data.frame()
colnames(df) = c("cell.types","samples","Freq")
df %<>% spread("samples","Freq")
rownames(df) = df$cell.types
df = df[order(df$cell.types),]

#write.csv(df, paste0(path,"Lung_24-",con,"_cell.types_by_samples.csv"))
write.csv(df, paste0(path,"Cell_types_by_samples.csv"))

# - Table: cluster name - cell type
df <- table(object$RNA_snn_res.0.8, object$cell.types) %>% 
        as.data.frame()
colnames(df) = c("RNA_snn_res.0.8","cell.types","Freq")
df %<>% spread("RNA_snn_res.0.8","Freq")
rownames(df) = df$cell.types
df = df[order(df$cell.types),]
write.csv(df, paste0(path,"Cluster_name-cell_type.csv"))

# - Table: number of cells per cluster (per each sample and total)
object$RNA_snn_res.0.8 %<>% as.character() %>% as.numeric()
df <- table(object$RNA_snn_res.0.8, object$orig.ident) %>% 
        as.data.frame()
colnames(df) = c("RNA_snn_res.0.8","samples","Freq")
df %<>% spread("samples","Freq")
rownames(df) = df$RNA_snn_res.0.8
df = df[order(as.numeric(df$RNA_snn_res.0.8)),]

#df = df[,-1]
#write.csv(df, paste0(path,"Lung_24-",con,"_cell.types_by_samples.csv"))
write.csv(df, paste0(path,"Clusters_by_samples.csv"))

#  UMAP, tSNE coordinates 
meta.data = cbind.data.frame(object@meta.data,
                             object@reductions$umap@cell.embeddings,
                             object@reductions$tsne@cell.embeddings)
meta.data = meta.data[,c("UMAP_1","UMAP_2","tSNE_1","tSNE_2","RNA_snn_res.0.8",
                         "cell.types","orig.ident","conditions")]
colnames(meta.data) %<>% sub("orig.ident","samples",.)
colnames(meta.data) %<>% sub("conditions","regions",.)

meta.data$RNA_snn_res.0.8 = as.numeric(as.character(meta.data$RNA_snn_res.0.8))

meta.data = meta.data[order(as.numeric(meta.data$RNA_snn_res.0.8)),]
print(colnames(meta.data))
write.csv(meta.data[,-c(3:4)], 
          paste0(path,"Lung_24_UMAP_coordinates.csv"))
write.csv(meta.data[,-c(1:2)], 
          paste0(path,"Lung_24_tSNE_coordinates.csv"))

# For P, D, T separately =============
Idents(object) = "conditions"

for (con in  c("proximal","distal","terminal")){
        sub_object <- subset(object,idents = con)
        #meta
        meta.data = cbind.data.frame(sub_object@meta.data,
                                     sub_object@reductions$umap@cell.embeddings,
                                     sub_object@reductions$tsne@cell.embeddings)
        meta.data = meta.data[,c("UMAP_1","UMAP_2","tSNE_1","tSNE_2","RNA_snn_res.0.8",
                                 "cell.types","orig.ident","conditions")]
        colnames(meta.data) %<>% sub("orig.ident","samples",.)
        colnames(meta.data) %<>% sub("conditions","regions",.)
        meta.data$RNA_snn_res.0.8 = as.numeric(as.character(meta.data$RNA_snn_res.0.8))
        
        meta.data = meta.data[order(as.numeric(meta.data$RNA_snn_res.0.8)),]
        print(colnames(meta.data))
        write.table(meta.data, paste0(path,"Lung_24-",con,"_meta.data.txt"), sep='\t', quote=F)
        
        meta_data <- cbind(rownames(sub_object@meta.data), sub_object@meta.data[,"cell.types", drop=F])   #####  cluster is the user’s specific cluster column
        write.table(meta_data, paste0(path,"Lung_24-",con,"_meta.data_cellphone.txt"), sep='\t', quote=F, row.names=F)
        
        # counts
        write.table(sub_object@assays$RNA@data, paste0(path,"Lung_24-",con,"_counts.txt"), sep='\t', quote=F)
        
}

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
df = df[order(as.numeric(df$cluster_res.0.8)),]
#write.csv(df,paste0(path,"Lung_24-",con,"_nGene_by_samples.csv"))
write.csv(df,paste0(path,"nGene_by_samples.csv"))

path <- "Yang/proximal_distal_terminal/Non-Integration/Counts and meta/"
terminal_meta_data <- read.csv(paste0(path,"Lung_24-proximal_meta.data.csv"),row.names = 1)
terminal_counts <- read.csv(paste0(path,"Lung_24-proximal_counts.csv"),row.names = 1)