########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","magrittr","MAST","data.table","biomaRt"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

save.path = "Yang/Lung_30/Cell_Phone_DB/"
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

#===== load data ===================
object = readRDS(file = "data/Lung_30_20200710.rds") 
print(unique(object@meta.data$conditions))
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
DefaultAssay(object) = "SCT"

#===== convert gene name to ENSEMBL===================


mart <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

ens <- rownames(object[["SCT"]]@data)
annotLookup <- getBM(mart=ensembl,
                     attributes=c("ensembl_gene_id","external_gene_name"),
                     filter="external_gene_name",
                     values=ens, uniqueRows=TRUE)
gene.keep = ens %in% annotLookup$external_gene_name
annotLookup = annotLookup[match(ens[gene.keep], annotLookup$external_gene_name),]


# For all
meta.data = cbind.data.frame(object@meta.data,
                             object@reductions$umap@cell.embeddings)
meta.data = meta.data[,c("UMAP_1","UMAP_2","SCT_snn_res.0.8",
                         "annotations3","orig.ident","conditions")]
colnames(meta.data) %<>% sub("orig.ident","samples",.)
colnames(meta.data) %<>% sub("conditions","regions",.)
colnames(meta.data) %<>% sub("annotations3","cell.types",.)
meta.data$SCT_snn_res.0.8 = as.numeric(as.character(meta.data$SCT_snn_res.0.8))
meta.data = meta.data[order(as.numeric(meta.data$SCT_snn_res.0.8)),]

write.csv(meta.data, paste0(save.path,"Cordinates_P_D_T_COPD.csv"),quote=F)

# For P, D, T separately =============
Idents(object) = "conditions"
#object$annotations3 %<>% gsub("-.*","",.) %>% gsub("[0-9]+","",.)
conditions <- c("proximal","distal","terminal","COPD")
for (con in  conditions){
        message(paste("write meta.data for", con))
        sub_object <- subset(object,idents = con)
        colnames(sub_object@meta.data) %<>% sub("annotations3","cell.types",.)
        
        #meta
        meta.data = cbind.data.frame(sub_object@meta.data,
                                     sub_object@reductions$umap@cell.embeddings)
        meta.data = meta.data[,c("UMAP_1","UMAP_2","SCT_snn_res.0.8",
                                 "cell.types","orig.ident","conditions")]
        colnames(meta.data) %<>% sub("orig.ident","samples",.)
        colnames(meta.data) %<>% sub("conditions","regions",.)

        meta.data$SCT_snn_res.0.8 = as.numeric(as.character(meta.data$SCT_snn_res.0.8))
        
        meta.data = meta.data[order(as.numeric(meta.data$SCT_snn_res.0.8)),]
        #write.table(meta.data, paste0(save.path,"Lung_30-",con,"_meta.data.txt"), sep='\t', quote=F)
        
        meta_data <- cbind("Cell" = rownames(sub_object@meta.data), sub_object@meta.data[,"cell.types", drop=F])   #####  cluster is the user’s specific cluster column
        write.table(meta_data, paste0(save.path,"Lung_30-",con,"_meta.data_cellphone.txt"), sep='\t', quote=F, row.names=F)
        
        message(paste("write counts for", con))
        counts = as.data.table(sub_object@assays$SCT@data,keep.rownames=TRUE)
        rownames(counts) = rownames(sub_object@assays$SCT@data)
        colnames(counts)[1] = "Gene"
        counts = counts[gene.keep,]
        print(table(counts$Gene == annotLookup$external_gene_name))
        counts$Gene = annotLookup$ensembl_gene_id
        
        fwrite(counts,
               file = paste0(save.path,"Lung_30-",con,"_counts.txt"),
               sep = "\t",row.names = FALSE, col.names = TRUE)
        Progress(which(conditions %in% con) , length(conditions))
        }

object_list <- SplitObject(object, split.by = "conditions")
counts_list <- lapply(object_list, function(x) x@assays$SCT@counts)
dt_list <- lapply(counts_list, function(x) as.data.table(x,keep.rownames=TRUE))
(conditions <- names(counts_list))

for(i in seq_along(conditions)){
        fwrite(dt_list[[i]],file = paste0(save.path,"Lung_30-",conditions[i],"_counts.txt"),
               sep = "\t",row.names = FALSE, col.names = TRUE)
        Progress(i, length(conditions))
}