########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("R.utils","Seurat","dplyr","kableExtra","ggplot2","scater",
        "scran","scran","BiocSingular"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
        }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
#args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0) stop("Sample info must be provided in excel!")
df_samples <- readxl::read_excel("doc/20200612_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",16))
#sample_n = intersect(sample_n, grep("COPD",df_samples$group,invert = T))
df_samples <- df_samples[sample_n,]
print(df_samples)
(samples = df_samples$sample)
nrow(df_samples)
# check missing data
current <- list.files("data")
(current <- current[!grepl(".Rda|RData",current)])
(missing_data <- df_samples$sample.id[!(df_samples$sample.id %in% current)])

# select species
if(unique(df_samples$species) %in% c("Homo_sapiens","Human")) species <- "hg19"
if(unique(df_samples$species) %in% c("Mus_musculus","Mouse")) species <- "mm10"
if(unique(df_samples$species) %in% c("Danio_rerio")) species <- "danRer10"

message("Copying the datasets")
if(length(missing_data)>0){
        # Move files from Download to ./data and rename them
        for(missing_dat in missing_data){
                old.pth  <- paste("~/Downloads", missing_dat,"outs",
                                  "filtered_feature_bc_matrix",sep = "/")
                list.of.files <- list.files(old.pth)
                new.folder <- paste("./data", missing_dat,"outs",
                                    "filtered_feature_bc_matrix",sep = "/")
                if(!dir.exists(new.folder)) dir.create(new.folder, recursive = T)
                # copy the files to the new folder
                file.copy(paste(old.pth, list.of.files, sep = "/"), new.folder)
                print(list_files <- list.files(new.folder))
        }
}
message("Loading the datasets")
Seurat_raw <- list()
Seurat_list <- list()

## Load the new version 10X dataset
sample_n = which(df_samples$tests %in% c("control",paste0("test",16)))
for(i in sample_n){
        Seurat_raw[[i]] <- Read10X(data.dir = paste0("data/",df_samples$sample.id[i],
                                                     "/outs/filtered_feature_bc_matrix/"))
        colnames(Seurat_raw[[i]]) = paste0(df_samples$sample[i],"_",colnames(Seurat_raw[[i]]))
        Seurat_list[[i]] <- CreateSeuratObject(Seurat_raw[[i]],
                                               min.cells = 0,
                                               min.features = 0)
        Seurat_list[[i]]@meta.data$tests <- df_samples$tests[i]
        
}
remove(Seurat_raw);GC()
#======1.1.2 record data quality before removing low quanlity cells =========================
# if args 2 is passed

message("QC")
cell.number <- sapply(Seurat_list, function(x) length(colnames(x)))
raw_nCount_RNA <- sapply(Seurat_list, function(x) mean(x$nCount_RNA))
raw_nFeature_RNA <- sapply(Seurat_list, function(x) mean(x$nFeature_RNA))

QC.list <- cbind(df_samples,cell.number, raw_nCount_RNA, raw_nFeature_RNA,
                 row.names = df_samples$sample)
write.csv(QC.list,paste0(path,"test31_QC_list.csv"))
QC.list %>% kable() %>% kable_styling()
remove(QC.list,cell.number,raw_nCount_RNA,raw_nFeature_RNA);GC()

#========1.1.3 g1 QC plots before filteration=================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), Seurat_list)
remove(Seurat_list);GC()


# read and select mitochondial genes
if(species == "hg19") mito = "^MT-"
if(species == "mm10") mito = "^mt-" # not Mt-
if(species == "danRer10") mito = "^mt-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)
Idents(object) = factor(Idents(object),levels = df_samples$sample)
g1 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
        VlnPlot(object = object, features = features, ncol = 3, pt.size = 0.01)+
                theme(axis.text.x = element_text(size=12),legend.position="none")
})
save(g1,file= paste0(path,"g1","_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))

#========1.2 scatter and scran ===============================
Seurat_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()

for(i in 1:length(df_samples$sample)){
        high.mito <- isOutlier(Seurat_list[[i]]$percent.mt, nmads=3, type="higher")
        low.lib <- isOutlier(log10(Seurat_list[[i]]$nCount_RNA), type="lower", nmad=3)
        low.genes <- isOutlier(log10(Seurat_list[[i]]$nFeature_RNA), type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        print(data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib), 
                   LowNgenes=sum(low.genes),Discard=sum(discard)))
        Seurat_list[[i]] <- Seurat_list[[i]][,!discard]
        #print(summary(!discard))
        print(i)
}

sce_list <- lapply(Seurat_list, as.SingleCellExperiment)
remove(Seurat_list);GC()

# cluster
set.seed(1000)
clusters_list <- lapply(sce_list,function(x){
        quickCluster(x, use.ranks=FALSE, BSPARAM=IrlbaParam())
}) 
sapply(clusters_list, table)

# computeSumFactors and normalize
sce_list <- mapply(function(x,y){
        computeSumFactors(x, min.mean=0.1, cluster=y)},
        x=sce_list,y=clusters_list)
lapply(sce_list,function(x) summary(sizeFactors(x)))
remove(clusters_list);GC()
#plot(sce_list[[1]]$nCount_RNA, sizeFactors(sce_list[[1]]), log="xy")
sce_list <- lapply(sce_list, normalize)

save(sce_list, file = paste0("data/","sce_",length(df_samples$sample),"_",gsub("-","",Sys.Date()),".Rda"))
