########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","magrittr","MAST"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#===== load data ===================
(load(file = "data/Lung_28_20200116.Rda"))
print(unique(object@meta.data$conditions))
# For P, D, T separately =============
Idents(object) = "conditions"
object$cell_types %<>% gsub("-.*","",.) %>% gsub("[0-9]+","",.)

for (con in  c("proximal","distal","terminal","COPD")){
        sub_object <- subset(object,idents = con)
        #meta
        meta.data = cbind.data.frame(sub_object@meta.data,
                                     sub_object@reductions$umap@cell.embeddings,
                                     sub_object@reductions$tsne@cell.embeddings)
        meta.data = meta.data[,c("UMAP_1","UMAP_2","tSNE_1","tSNE_2","SCT_snn_res.0.8",
                                 "cell_types","orig.ident","conditions")]
        colnames(meta.data) %<>% sub("orig.ident","samples",.)
        colnames(meta.data) %<>% sub("conditions","regions",.)
        meta.data$SCT_snn_res.0.8 = as.numeric(as.character(meta.data$SCT_snn_res.0.8))
        
        meta.data = meta.data[order(as.numeric(meta.data$SCT_snn_res.0.8)),]
        print(colnames(meta.data))
        write.table(meta.data, paste0(path,"Lung_28-",con,"_meta.data.txt"), sep='\t', quote=F)
        
        meta_data <- cbind(rownames(sub_object@meta.data), sub_object@meta.data[,"cell.types", drop=F])   #####  cluster is the user’s specific cluster column
        write.table(meta_data, paste0(path,"Lung_28-",con,"_meta.data_cellphone.txt"), sep='\t', quote=F, row.names=F)
}

# - Table: number of expressed genes per cluster (per each sample and total) ==============
df <- table(object$SCT_snn_res.0.8, object$orig.ident) %>% as.data.frame()
colnames(df) =c("cluster_res.0.8","samples","nGene")

df$cluster_res.0.8 <- as.character(df$cluster_res.0.8)
df$samples <- as.character(df$samples)

for(i in seq_len(nrow(df))) {
        cells <- object$SCT_snn_res.0.8 %in% df[i,"cluster_res.0.8"] & object$orig.ident %in% df[i,"samples"]
        df[i,"nGene"] = as.integer(mean(object@meta.data[cells,"nFeature_RNA"]))
        svMisc::progress(i/nrow(df)*100)
}
df %<>% spread("samples","nGene")
df[is.na(df)] = 0
df = df[order(as.numeric(df$cluster_res.0.8)),]
#write.csv(df,paste0(path,"Lung_28-",con,"_nGene_by_samples.csv"))
write.csv(df,paste0(path,"nGene_by_samples.csv"))