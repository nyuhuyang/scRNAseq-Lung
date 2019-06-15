########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(DropletUtils)
library(dplyr)
library(kableExtra)
library(scater)
#install_github("MarioniLab/scran") #BiocManager::install("BiocNeighbors", version = "devel")
library(scran)
library(EnsDb.Hsapiens.v86)
library(devtools)
library(Matrix)
library(devtools)
#library(scRNAseq)#BiocInstaller::biocLite("scRNAseq")
source("../R/Seurat_functions.R")
source("../R/scatter_utils.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
########################################################################
#
#  0.1~0.5 scater
# 
# ######################################################################
# 0.1. Setting up the data
# read sample summary list
df_samples <- readxl::read_excel("doc/20190119_Single_cell_sample list.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c(paste0("test",4)))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
sample.id <- df_samples$sample.id[sample_n]
conditions <- df_samples$conditions[sample_n]
projects <- df_samples$project[sample_n]
tests <- df_samples$tests[sample_n] 

sce_list <- list()
species <- "hg19"
for(i in 1:length(samples)){
        fname <- paste0("./data/",sample.id[i],
                        "/outs/filtered_feature_bc_matrix")
        sce_list[[i]] <- read10xCounts.1(fname, col.names=TRUE,
                                         add.colnames = samples[i])
}
names(sce_list) <- samples
# 0.1.2 Annotating the rows
for(i in 1:length(samples)){
        rownames(sce_list[[i]]) <- uniquifyFeatureNames(rowData(sce_list[[i]])$ID,
                                                        rowData(sce_list[[i]])$Symbol)
        print(head(rownames(sce_list[[i]]),3))
        print(length(rownames(sce_list[[i]])))
}

# We also identify the chromosomal location for each gene. 
# The mitochondrial percentage is particularly useful for later quality control.
for(i in 1:length(samples)){
        location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce_list[[i]])$ID, 
                           column="SEQNAME", keytype="GENEID")
        rowData(sce_list[[i]])$CHR <- location
        print(summary(location=="MT"))
}

# 0.4 Quality control on the cells#########################
# It is entirely possible for droplets to contain damaged or dying cells,
# which need to be removed prior to downstream analysis. 
# We compute some QC metrics using  calculateQCMetrics() (McCarthy et al. 2017) 
# and examine their distributions in Figure 2.
sce_list <- lapply(sce_list, function(x) calculateQCMetrics(x,compact = FALSE,
                        feature_controls=list(Mito=which(location=="MT"))))

########################################################################

# Ideally, we would remove cells with low library sizes or total number of expressed features as described previously.
# However, this would likely remove cell types with low RNA content,
# especially in a heterogeneous population with many different cell types.
# Thus, we use a more relaxed strategy and only remove cells with large mitochondrial proportions,
# using it as a proxy for cell damage. 
# (Keep in mind that droplet-based datasets usually do not have spike-in RNA.)
# Low-quality cells are defined as those with extreme values for these QC metrics and are removed.
for(i in 1:length(samples)){
        high.mito <- isOutlier(sce_list[[i]]$pct_counts_Mito, nmads=3, type="higher")
        low.lib <- isOutlier(sce_list[[i]]$log10_total_counts, type="lower", nmad=3)
        low.genes <- isOutlier(sce_list[[i]]$log10_total_features_by_counts, type="lower", nmad=3)
        discard <- high.mito | low.lib | low.genes
        data.frame(HighMito= sum(high.mito),LowLib=sum(low.lib), 
                   LowNgenes=sum(low.genes),Discard=sum(discard))
        sce_list[[i]] <- sce_list[[i]][,!discard]
        print(summary(!discard))
}

########################################################################
#
#  0.6~ scran
# 
# ######################################################################
# Use natural Log transform to fit Seurat

for(i in 1:length(sce_list)){
        logcounts(sce_list[[i]]) <- as(log1p(assay(sce_list[[i]], "counts")),"dgCMatrix")
}
# 0.6 Normalizing for cell-specific biases
#clusters <- list()
#for(i in 1:length(sce_list)){
#        clusters[[i]] <- quickCluster(sce_list[[i]], method="igraph", min.size=50,
#                                      assay.type = "logcounts",
#                                 irlba.args=list(maxit=1000)) # for convergence.
#        print(table(clusters[[i]]))
#        sce_list[[i]] <- computeSumFactors(sce_list[[i]], min.mean=0.1, 
#                                           cluster=clusters[[i]])
#        print(summary(sizeFactors(sce_list[[i]])))
#        print(paste0(i,":",length(sce_list)," done"))
#}

save(sce_list, file = "data/sce_5_20190123.Rda")
