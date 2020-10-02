# How to download dbGaP CF Genotype-Tissue Expression Project (GTEx) data set
# https://gtexportal.org/home/datasets
# Download https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz"
# Download https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt"

devtools::install_github('xuranw/MuSiC')
library(MuSiC)
library(Biobase)
library(xbioc)
library(reshape2)
library(cowplot)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ======= load seurat object ================
(load(file="data/Lung_GTEx_20200307.Rda"))
bulk.eset <- new("ExpressionSet", exprs = as.matrix(object@assays$SCT@data),
                 phenoData = AnnotatedDataFrame(object@meta.data))

object = readRDS(file="data/Lung_28_Global_20200219.rds")

object_list <- SplitObject(object, split.by = "conditions")

sc.eset_list <- list()
for(i in seq_along(object_list)){
        sc.eset_list[[i]] = new("ExpressionSet", exprs = as.matrix(object_list[[i]]@assays$SCT@data),
                                 phenoData = AnnotatedDataFrame(object_list[[i]]@meta.data))
        Progress(i, length(object_list))
}

sc.eset <- sc.eset_list[[1]]
for(i in 2:length(sc.eset_list)){
        sc.eset <- BiocGenerics::combine(sc.eset, sc.eset_list[[i]])
        Progress(i, length(sc.eset_list))
}
ncol(object) == ncol(sc.eset)
remove(object_list,sc.eset_list); GC()

Idents(object) = "cell.labels"
Epi <- subset(object, idents = c("AT-p","AT1","AT2","BC","BC-p","C1",
                                 "C2","C3","C4","Cr",
                                 "H","IC","Ion","MEC","NEC","Nr",
                                 "p-C", "Pr","S","S-d","SM1","SM2","SM3",
                                 "SMG-Muc","SMG-Ser","Sq"))
Epi.eset = new("ExpressionSet", exprs = as.matrix(Epi@assays$SCT@data),
                        phenoData = AnnotatedDataFrame(Epi@meta.data))
rm(Epi,object);GC()
# Estimate all cell type proportions
sc.eset$major_cell.labels = sc.eset$cell.labels
sc.eset$major_cell.labels %<>% gsub("[1-9]$","",.)
sc.eset$major_cell.labels %<>% gsub("^p-C","C",.)
sc.eset$major_cell.labels %<>% gsub("-.*","",.)
choose = c("all","major","no-Ion","no-Ion_major")[4]
if(choose == "all") cell.labels = sort(unique(sc.eset$cell.labels))
if(choose == "major") cell.labels = sort(unique(sc.eset$major_cell.labels))
if(choose == "no-Ion") cell.labels = sort(unique(sc.eset$cell.labels)) %>% .[-grep("Ion",.)]
if(choose == "no-Ion_major") cell.labels = sort(unique(sc.eset$major_cell.labels)) %>% .[-grep("Ion",.)]

Lung.prop = music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset,
                               clusters = ifelse(choose %in% c("all","no-Ion"),
                                                 "cell.labels","major_cell.labels"),
                               samples = 'orig.ident', 
                               select.ct = cell.labels, verbose = T)
str(Lung.prop)
names(Lung.prop)

# Jitter plot of estimated all cell type proportions
cell.labels = sort(colnames(Lung.prop$Est.prop.weighted))
Lung.prop$Est.prop.weighted = Lung.prop$Est.prop.weighted[,cell.labels]
jitter.fig = Jitter_Est(list(data.matrix(Lung.prop$Est.prop.weighted),
                             data.matrix(Lung.prop$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jpeg(paste0(path,"All_",choose,"_celltype.jpeg"), units="in", width=10, height=7,res=600)
print(jitter.fig)
dev.off()

# Estimate all cell type proportions
Epi.eset$major_cell.labels = Epi.eset$cell.labels
Epi.eset$major_cell.labels %<>% gsub("^C[1-4]","C",.)
Epi.eset$major_cell.labels %<>% gsub("^SM[1-4]","SM",.)
Epi.eset$major_cell.labels %<>% gsub("^p-C","C",.)
Epi.eset$major_cell.labels %<>% gsub("-.*","",.)
Epi.eset$major_cell.labels %<>% gsub("^AT$","AT2",.)
major_cell.labels = sort(unique(Epi.eset$major_cell.labels)) %>% .[-grep("Ion",.)]
.[-grep("Ion",major_cell.labels)]

choose = c("all","major","no-Ion","no-Ion_major")[4]
if(choose == "all") Epi.labels = sort(unique(Epi.eset$cell.labels))
if(choose == "major") Epi.labels = sort(unique(Epi.eset$major_cell.labels))
if(choose == "no-Ion") Epi.labels = sort(unique(Epi.eset$cell.labels)) %>% .[-grep("Ion",.)]
if(choose == "no-Ion_major") Epi.labels = sort(unique(Epi.eset$major_cell.labels)) %>% .[-grep("Ion",.)]

Epi.prop = music_prop(bulk.eset = bulk.eset, sc.eset = Epi.eset,
                       clusters = ifelse(choose %in% c("all","no-Ion"),
                                         "cell.labels","major_cell.labels"),
                       samples = 'orig.ident', 
                       select.ct = Epi.labels, verbose = T)
names(Epi.prop)

# Jitter plot of estimated Epi cell type proportions
jitter.fig = Jitter_Est(list(data.matrix(Epi.prop$Est.prop.weighted),
                             data.matrix(Epi.prop$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jpeg(paste0(path,"Epi_",choose,"_celltype.jpeg"), units="in", width=10, height=7,res=600)
print(jitter.fig)
dev.off()

# Estimation of cell type proportions #####################
# Clustering single cell data
# Produce the first step information
sc.eset

sc.basis = music_basis(sc.eset, clusters = 'cell.labels', samples = 'orig.ident', 
                             select.ct = sort(unique(sc.eset$cell.labels)))

d <- dist(t(log(sc.basis$Disgn.mtx + 1e-6)), method = "euclidean")
hc1 <- hclust(d, method = "ward.D" )
d <- dist(t(log(sc.basis$M.theta + 1e-8)), method = "euclidean")
hc2 <- hclust(d, method = "ward.D")

jpeg(paste0(path,"Celltypes_clusters.jpeg"), units="in", width=10, height=7,res=600)
par(mfrow = c(1, 2))
plot(hc1, cex = 0.5, hang = -1, main = 'Cluster log(Design Matrix)')
plot(hc2, cex = 0.5, hang = -1, main = 'Cluster log(Mean of RA)')
dev.off()

# Produce the first step information ====================
Epi.basis = music_basis(Epi.eset, clusters = 'RNA_snn_res.0.3', samples = 'orig.ident', 
                       select.ct = sort(unique(Epi.eset$RNA_snn_res.0.3)))

d <- dist(t(log(Epi.basis$Disgn.mtx + 1e-6)), method = "euclidean")
hc1 <- hclust(d, method = "complete" )
d <- dist(t(log(Epi.basis$M.theta + 1e-8)), method = "euclidean")
hc2 <- hclust(d, method = "complete")

jpeg(paste0(path,"Epi_cell_clusters~.jpeg"), units="in", width=10, height=7,res=600)
par(mfrow = c(1, 2))
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')
dev.off()

hc1_clust <- cutree(hc1, k = 5);sort(hc1_clust)
hc2_clust <- cutree(hc2, k = 5);sort(hc2_clust)
####################################
# Bulk Tissue Cell Type Estimation with Pre-grouping of Cell Types
###################################
clusters.type = list(cluster1 = names(hc1_clust)[hc1_clust == 1],
                     cluster2 = names(hc1_clust)[hc1_clust == 2],
                     cluster3 = names(hc1_clust)[hc1_clust == 3],
                     cluster4 = names(hc1_clust)[hc1_clust == 4],
                     cluster5 = names(hc1_clust)[hc1_clust == 5])

cl.type = as.character(sc.eset$cell.labels)

for(cl in 1:length(clusters.type)){
        cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
pData(sc.eset)$clusterType = factor(cl.type, levels = names(clusters.type))

# We manually specify the cluster and annotated single cell data with cluster information. 
table(pData(EMTAB.eset)$cellType)
(celltype = unique(pData(EMTAB.eset)$cellType) %>% as.character)

Est.mouse.bulk = music_prop.cluster(bulk.eset = Mouse.bulk.eset, 
                                    sc.eset = Mousesub.eset, 
                                    group.markers = IEmarkers, clusters = 'cellType', group = 'clusterType', samples = 'sampleID', clusters.type = clusters.type)

