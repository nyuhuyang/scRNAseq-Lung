library(SoupX)
library(magrittr)
library(Seurat) # Seurat 4
library(sctransform)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
#https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/202108014_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()

df_samples = df_samples[grepl("P-norm|D-norm|T-norm|D-COPD", df_samples$condition),]
df_samples = df_samples[!grepl("UNC_44_P|VU_29_D|VU_35_D", df_samples$sample),]

# check missing data
read.path = "data/scRNA-seq/counts"
current <- list.files(read.path)
(missing_data <- df_samples$sample.id[!(df_samples$sample.id %in% current)])

#============== filtered counts & re-clustering ====================
adj.matrix_list <- pbapply::pblapply(df_samples$sample.id, function(s){
    readDir <- file.path(read.path,as.character(s),"outs")
    filt.matrix <- Seurat::Read10X(file.path(readDir, "filtered_feature_bc_matrix"))
    raw.matrix <- Seurat::Read10X(file.path(readDir, "raw_feature_bc_matrix"))
    soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
    
    obj <- CreateSeuratObject(filt.matrix,min.cells = 0,names.delim = "-",min.features = 0) %>%
        NormalizeData(verbose = FALSE) %>%
        FindVariableFeatures(verbose = FALSE) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(verbose = FALSE) %>%
        FindNeighbors(reduction = "pca",dims = 1:30) %>%
        FindClusters(resolution = 0.8, algorithm = 1,verbose = F)
    
    
    soup.channel  <- setClusters(soup.channel, setNames(obj$RNA_snn_res.0.8, colnames(obj)))
    soup.channel  <- autoEstCont(soup.channel, priorRhoStdDev = 0.3)
    adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
    return(adj.matrix)
})

names(adj.matrix_list) = df_samples$sample
    
for (s in df_samples$sample) {
    colnames(adj.matrix_list[[s]]) %<>% gsub("-[0-9+]","",.)
    colnames(adj.matrix_list[[s]]) = paste0(s,"-", colnames(adj.matrix_list[[s]]))
}

adj.matrix <- do.call(cbind, adj.matrix_list)
meta.data = readRDS("output/Lung_30_20210831_metadata_v2.rds")
table(colnames(adj.matrix) %in% rownames(meta.data))
table(rownames(meta.data) %in% colnames(adj.matrix))

object <- CreateSeuratObject(adj.matrix[,rownames(meta.data)],
                             min.cells = 0,names.delim = ".",min.features = 0)
mito <- "^MT-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)

format(object.size(object)*60,unit = "GB")
options(future.globals.maxSize= object.size(object)*60)
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
object[["SCT"]]@scale.data = matrix(0,0,0)

object_orig = readRDS(file = "data/Lung_SCT_30_20210831.rds")
object[["umap"]]  <- object_orig[["umap"]]
object[["tsne"]]  <- object_orig[["tsne"]]

if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}

format(object.size(object),unit = "GB")
saveRDS(object, file = "data/Lung_30_SoupX_20221101.rds")
