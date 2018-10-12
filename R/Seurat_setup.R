########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(Matrix)
library(sva)
library(SingleR)
source("../R/Seurat_functions.R")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# Load the dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

# rename all "_" into "-" in sample names 
EC_raw <- list()
EC_Seurat <- list()
samples <- c("Black6-Control-1","Black6-Contorl-2")
conditions <- c("Control","Control")
for(i in 1:length(samples)){
        EC_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                   samples[i],"/outs/filtered_gene_bc_matrices/mm10/"))
        colnames(EC_raw[[i]]) <- paste0(samples[i],"_",colnames(EC_raw[[i]]))
        EC_Seurat[[i]] <- CreateSeuratObject(EC_raw[[i]],
                                               min.cells = 3,
                                               min.genes = 200,
                                               names.delim = "_",
                                               project = "paula")
        EC_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
EC <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), EC_Seurat)
remove(EC_raw,EC_Seurat);GC()
EC <- FilterCells(EC, subset.names = "nGene",
                    low.thresholds = 200,
                    high.thresholds = Inf) %>%
        NormalizeData() %>%
        ScaleData(display.progress = FALSE) %>%
        FindVariableGenes(do.plot = FALSE, display.progress = FALSE)
#save(EC, file = "./data/EC_20180825.Rda")

#======1.2 QC, pre-processing and normalizing the data=========================
# 1.2.1 Calculate median UMI per cell
EC_raw_data <- as.matrix(x = EC@raw.data)
mean(colSums(EC_raw_data))
median(colSums(EC_raw_data))
min(colSums(EC_raw_data))
remove(EC_raw_data);GC()

# 1.2.2 Identifying Outlier Cells
# Plot genes per cell
# How many genes expressed per cells
complexity.per.cell <- apply(EC@data, 2, function(x) sum(x>0))
# Mean count per cell.
mean.count.per.cell <- apply(EC@data, 2, function(x) mean(x))
# Gene prevalence
gene.prevalence <- apply(EC@raw.data, 1, function(x) sum(x>0))
# Complexity by mean expression
par(mfrow = c(2,1))
plot(complexity.per.cell, mean.count.per.cell)

# 1.2.3 calculate mitochondria percentage
mito.genes <- grep(pattern = "^mt-", x = rownames(x = EC@data), value = TRUE)
percent.mito <- Matrix::colSums(EC@raw.data[mito.genes, ])/Matrix::colSums(EC@raw.data)
EC <- AddMetaData(object = EC, metadata = percent.mito, col.name = "percent.mito")

g1 <- VlnPlot(object = EC, features.plot = c("percent.mito"), nCol = 1,
             x.lab.rot = T, do.return = T)
g1
EC <- FilterCells(object = EC, subset.names = c("nUMI"),
                    low.thresholds = c(70000), 
                    high.thresholds = c(Inf))

g2 <- VlnPlot(object = EC, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 1,
              x.lab.rot = T, do.return = T)

par(mfrow = c(2, 1))
GenePlot(object = EC, gene1 = "nUMI", gene2 = "percent.mito",use.raw = T)
GenePlot(object = EC, gene1 = "nUMI", gene2 = "nGene",use.raw = T)

VlnPlot(object = EC, features.plot = c("nUMI"), nCol = 1,
              group.by = "ident", x.lab.rot = T, do.return = T)
plot_grid(g1,g2)
# After removing unwanted cells from the dataset, the next step is to normalize the data.
EC <- NormalizeData(object = EC, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
EC <- FindVariableGenes(object = EC, mean.function = ExpMean, 
                          dispersion.function = LogVMR, do.plot = FALSE, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(EC@var.genes)
#======1.3 1st run of pca-tsne  =========================
EC <- ScaleData(object = EC) %>%
        RunPCA() %>%
        FindClusters(dims.use = 1:20, force.recalc = T, print.output = FALSE) %>%
        RunTSNE()
#EC@meta.data$orig.ident <- gsub("PND18pre","PND18",EC@meta.data$orig.ident)
TSNEPlot(object = EC, do.label = F, group.by = "orig.ident", 
         do.return = TRUE, no.legend = F,# colors.use = singler.colors,
         pt.size = 1,label.size = 8 )+
        ggtitle("TSNEPlot of all samples")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

save(EC, file = "./data/EC_20180825.Rda")
Iname = load("./data/EC_20180825.Rda")

#======1.4 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "./data/seurat_resources/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- MouseGenes(EC,cc.genes[1:43])
g2m.genes <- MouseGenes(EC,cc.genes[44:97])
# Assign Cell-Cycle Scores
EC <- CellCycleScoring(object = EC, s.genes = s.genes, g2m.genes = g2m.genes, 
                         set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = EC, features.plot = MouseGenes(EC,c("PCNA", "TOP2A", "MCM6", "MKI67")), 
          nCol = 2)
# regressing out the difference between the G2M and S phase scores
EC@meta.data$CC.Difference <- EC@meta.data$S.Score - EC@meta.data$G2M.Score
# view cell cycle scores and phase assignments
head(x = EC@meta.data)

#======1.5 Add batch id =========================
batchname = EC@meta.data$orig.ident
batch.effect = rep(NA,length(batchname))
batch.effect[batchname %in% c("PND14","PND18","PND18pre","PND30")] = 1
batch.effect[batchname %in% c("PND06","PND25","Ad-depleteSp","Ad-Thy1")] = 2
names(batch.effect) = rownames(EC@meta.data)
EC <- AddMetaData(object = EC, metadata = batch.effect, col.name = "batch.effect")
table(EC@meta.data$batch.effect)
head(x = EC@meta.data)
#======1.6 batch-correct using ComBat =========================
SingleFeaturePlot.1(EC,"nUMI",threshold=15000)
SingleFeaturePlot.1(EC,"batch.effect",threshold=1.0)
SingleFeaturePlot.1(EC,"percent.mito",threshold=0.05)
SingleFeaturePlot.1(EC,"CC.Difference",threshold=0.05)
m = as.matrix(EC@data)
m = m[rowSums(m)>0,]
com = ComBat(m, batch.effect, prior.plots=FALSE, par.prior=TRUE)
#----save files just in case------
saveRDS(Matrix(as.matrix(com)), file = "./data/Combat_data.Rda")
saveRDS(EC@data, file = "./data/EC_data.Rda")
remove(m,com);GC()
#---------------------
EC@data = readRDS("./data/Combat_data.Rda")
EC@scale.data = NULL
EC <- ScaleData(object = EC,#genes.use = EC@var.genes,
                  model.use = "negbinom", do.par=T, do.center = T, do.scale = T,
                  vars.to.regress = c("CC.Difference"),#"CC.Difference","percent.mito"--nogood,"nUMI"--nogood
                  display.progress = T)
EC@data =  readRDS("./data/EC_data.Rda")

gene.use <- rownames(EC@scale.data);length(gene.use)
EC@data = EC@data[gene.use,]
#save(EC, file = "./data/EC_20181001.Rda") #do.center = F, do.scale = T
#======1.7 unsupervised clustering =========================
EC <- RunPCA(object = EC, pc.genes = EC@var.genes, pcs.compute = 100, 
               do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = EC)
PCElbowPlot(object = EC, num.pc = 100)
PCHeatmap(EC, pc.use = c(1:3, 25:30), cells.use = 500, do.balanced = TRUE)
EC <- RunTSNE(object = EC, reduction.use = "pca", dims.use = 1:30, 
                do.fast = TRUE, perplexity= 30)

EC <- FindClusters(object = EC, reduction.type = "pca", dims.use = 1:30, 
                     resolution = 0.8, 
                     k.param = 30,force.recalc = T,
                     save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)
#EC@meta.data$orig.ident <- gsub("PND18pre","PND18",EC@meta.data$orig.ident)
TSNEPlot.1(object = EC, do.label = T, group.by = "ident", 
         do.return = TRUE, no.legend = T, 
         text.repel = T, label.repel = F,
        #colors.use = singler.colors,
         pt.size = 1,label.size = 6 )+
        ggtitle("TSNEPlot, resolution = 0.8, k.param = 30")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 
EC <- StashIdent(object = EC, save.name = "res.0.8")
save(EC, file = "./data/EC_20180825.Rda")
