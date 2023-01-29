library(SoupX)
library(magrittr)
library(Seurat) # Seurat 4
library(sctransform)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
#https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html
#======1.1 Setup the Seurat objects =========================
# read sample summary list
# read sample summary list
df_samples <- readxl::read_excel("output/20220523/20220523_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
nrow(df_samples)
df_samples$date %<>% gsub(" UTC","",.) %>% as.character()
#============== filtered counts & re-clustering ====================
adj.matrix_list <- pbapply::pblapply(df_samples$sample.id, function(s){
    readDir <- file.path(read.path,as.character(s),"outs")
    filt.matrix <- Seurat::Read10X(file.path(readDir, "filtered_feature_bc_matrix"),strip.suffix = TRUE)
    raw.matrix <- Seurat::Read10X(file.path(readDir, "raw_feature_bc_matrix"),strip.suffix = TRUE)
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
    colnames(adj.matrix_list[[s]]) = paste0(s,"-", colnames(adj.matrix_list[[s]]))
}


adj.matrix <- do.call(cbind, adj.matrix_list)
meta.data <- readRDS("output/Lung_time15_metadata_20220523_v2.rds")
table(colnames(adj.matrix) %in% rownames(meta.data))
table(rownames(meta.data) %in% colnames(adj.matrix))

object <- CreateSeuratObject(adj.matrix[,rownames(meta.data)],
                             min.cells = 0,names.delim = ".",min.features = 0)
mito <- "^MT-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)


object_orig = readRDS(file = "data/Lung_time15_20220523.rds")
object@commands <- object_orig@commands

format(object.size(object)*15,unit = "GB")
options(future.globals.maxSize= object.size(object)*15)
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
object[["SCT"]]@scale.data = matrix(0,0,0)

object@reductions  <- object_orig@reductions

if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}

format(object.size(object),unit = "GB")
saveRDS(object, file = "data/Lung_time15_SoupX_20230129.rds")
