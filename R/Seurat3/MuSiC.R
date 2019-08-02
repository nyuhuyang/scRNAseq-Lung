# How to download dbGaP CF Genotype-Tissue Expression Project (GTEx) data set
# https://gtexportal.org/home/datasets
# Download https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz"
# Download https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt"

#devtools::install_github('xuranw/MuSiC')
library(MuSiC)
library(Biobase)
library(xbioc)
library(reshape2)
library(cowplot)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ======= load seurat object ================
(load(file="output/Lung_bulk_427_20190731.Rda"))
bulk.eset <- new("ExpressionSet", exprs = as.matrix(object@assays$SCT@data),
                 phenoData = AnnotatedDataFrame(object@meta.data))

(load(file="data/Lung_harmony_12_20190614.Rda"))
sc.eset <- new("ExpressionSet", exprs = as.matrix(object@assays$RNA@data),
                 phenoData = AnnotatedDataFrame(object@meta.data))

(load("data/Epi_harmony_12_20190627.Rda"))
Epi.eset <- new("ExpressionSet", exprs = as.matrix(Epi@assays$RNA@data),
               phenoData = AnnotatedDataFrame(Epi@meta.data))

remove(object,Epi); GC()
# Estimate all cell type proportions
Lung.prop = music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset,
                               clusters = 'cell.type',
                               samples = 'orig.ident', 
                               select.ct = unique(sc.eset$cell.type), verbose = T)
str(Lung.prop)
names(Lung.prop)

# Jitter plot of estimated all cell type proportions
cell.type = sort(colnames(Lung.prop$Est.prop.weighted))
Lung.prop$Est.prop.weighted = Lung.prop$Est.prop.weighted[,cell.type]
jitter.fig = Jitter_Est(list(data.matrix(Lung.prop$Est.prop.weighted),
                             data.matrix(Lung.prop$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jpeg(paste0(path,"jitter_fig_celltype_prop.jpeg"), units="in", width=10, height=7,res=600)
print(jitter.fig)
dev.off()

# Estimate all cell type proportions
Epi.prop = music_prop(bulk.eset = bulk.eset, sc.eset = Epi.eset,
                       clusters = 'cell.type',
                       samples = 'orig.ident', 
                       select.ct = unique(Epi.eset$cell.type), verbose = T)
names(Epi.prop)

# Jitter plot of estimated Epi cell type proportions
jitter.fig = Jitter_Est(list(data.matrix(Epi.prop$Est.prop.weighted),
                             data.matrix(Epi.prop$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jpeg(paste0(path,"jitter_fig_Epi_celltype_prop.jpeg"), units="in", width=10, height=7,res=600)
print(jitter.fig)
dev.off()

# Produce the first step information ====================
sc.basis = music_basis(sc.eset, clusters = 'cell.type', samples = 'orig.ident', 
                             select.ct = sort(unique(sc.eset$cell.type)))

d <- dist(t(log(sc.basis$Disgn.mtx + 1e-6)), method = "euclidean")
hc1 <- hclust(d, method = "complete" )
d <- dist(t(log(sc.basis$M.theta + 1e-8)), method = "euclidean")
hc2 <- hclust(d, method = "complete")

jpeg(paste0(path,"All_cell_clusters.jpeg"), units="in", width=10, height=7,res=600)
par(mfrow = c(1, 2))
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')
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

####################################
# Bulk Tissue Cell Type Estimation with Pre-grouping of Cell Types
###################################
clusters.type = list(C1 = 'Neutro', C2 = 'Podo', C3 = c('Endo', 'CD-PC', 'LOH', 'CD-IC', 'DCT', 'PT'), C4 = c('Macro', 'Fib', 'B lymph', 'NK', 'T lymph'))

cl.type = as.character(Mousesub.eset$cellType)

for(cl in 1:length(clusters.type)){
        cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
pData(Mousesub.eset)$clusterType = factor(cl.type, levels = c(names(clusters.type), 'CD-Trans', 'Novel1', 'Novel2'))


# We manually specify the cluster and annotated single cell data with cluster information. 
table(pData(EMTAB.eset)$cellType)
(celltype = unique(pData(EMTAB.eset)$cellType) %>% as.character)

Est.mouse.bulk = music_prop.cluster(bulk.eset = Mouse.bulk.eset, 
                                    sc.eset = Mousesub.eset, 
                                    group.markers = IEmarkers, 
                                    clusters = 'cellType', 
                                    group = 'clusterType', 
                                    samples = 'sampleID', 
                                    clusters.type = celltype)

