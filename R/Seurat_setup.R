########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(scran)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/181012_Single_cell_sample list.xlsx")
samples <- df_samples$samples
projects <- df_samples$projects
conditions <- df_samples$conditions
df_samples %>% kable() %>% kable_styling()
samples_ddSeq <- df_samples$samples[df_samples$projects %in% "ddSeq"]
samples_10X <- df_samples$samples[df_samples$projects %in% "EC-SR-5444"]

Lung_raw <- list()
Lung_Seurat <- list()
for(i in 1:length(samples_ddSeq)){
    Lung_raw[[i]] <- read.table(file = paste0("./data/",samples_ddSeq[i],
                                              "/counts.merged.txt.gz"),
                                header = T, row.names =1)
    colnames(Lung_raw[[i]]) <- paste0(samples_ddSeq[i],
                                      "_",colnames(Lung_raw[[i]]))
}

for(i in (length(samples_ddSeq)+1):length(samples)){
    Lung_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                             samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
    colnames(Lung_raw[[i]]) <- paste0(samples[i],
                                            "_",colnames(Lung_raw[[i]]))
}

for(i in 1:length(samples)){
    Lung_Seurat[[i]] <- CreateSeuratObject(Lung_raw[[i]],
                                                 min.cells = 3,
                                                 min.genes = 0,
                                                 project = projects[i],
                                                 names.delim = "_")
    Lung_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
#======1.1.2 QC before merge =========================
cell.number <- sapply(Lung_Seurat, function(x) length(x@cell.names))
QC_list <- lapply(Lung_Seurat, function(x) as.matrix(x = x@raw.data))
median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))

min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))

QC.list <- cbind(df_samples,cell.number, median.nUMI,median.nGene,min.nUMI,min.nGene,
                      row.names = samples)
write.csv(QC.list,paste0(path,"QC_list.csv"))
QC.list %>% kable() %>% kable_styling()
remove(QC_list);GC()
#========1.1.3 merge ===================================
Lung <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), Lung_Seurat)
remove(Lung_raw,Lung_Seurat);GC()
Lung <- FilterCells(Lung, subset.names = "nGene",
                    low.thresholds = 50,
                    high.thresholds = Inf) %>%
    NormalizeData() %>%
    ScaleData(display.progress = FALSE) %>%
    FindVariableGenes(do.plot = FALSE, display.progress = FALSE)
save(Lung, file = "./data/Lung_20181101.Rda")

#======1.2 QC, pre-processing and normalizing the data=========================
# 1.2.1 Calculate median UMI per cell
Iname = load(file = "./data/Lung_20181101.Rda")
# 1.2.3 calculate mitochondria percentage
mito.genes <- grep(pattern = "^MT-", x = rownames(x = Lung@data), value = TRUE)
percent.mito <- Matrix::colSums(Lung@raw.data[mito.genes, ])/Matrix::colSums(Lung@raw.data)
Lung <- AddMetaData(object = Lung, metadata = percent.mito, col.name = "percent.mito")
#Lung@ident = factor(Lung@ident,levels = samples)

g1 <- VlnPlot(object = Lung, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

Lung <- FilterCells(object = Lung, subset.names = c("nGene","nUMI","percent.mito"),
                    low.thresholds = c(200,300, -Inf), 
                    high.thresholds = c(Inf,Inf,0.1))

g2 <- VlnPlot(object = Lung, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

jpeg(paste0(path,"/S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                    scale_y_log10(limits = c(100,7000)),#+ylim(c(0,1000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                    scale_y_log10(limits = c(100,7000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                    scale_y_log10(limits = c(200,60000)),#+ylim(c(0,1000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                    scale_y_log10(limits = c(200,60000))))
dev.off()
jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                    ylim(c(0,0.1)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,0.1))))
dev.off()
######################################

# After removing unwanted cells from the dataset, the next step is to normalize the data.
Lung <- NormalizeData(object = Lung, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
Lung <- FindVariableGenes(object = Lung, mean.function = ExpMean, 
                          dispersion.function = LogVMR, do.plot = FALSE, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(Lung@var.genes)
#======1.3 1st run of pca-tsne  =========================
Lung <- ScaleData(object = Lung) %>%
    RunPCA() %>%
    FindClusters(dims.use = 1:20, force.recalc = T, print.output = FALSE) %>%
    RunTSNE()

p1 <- TSNEPlot(object = Lung, do.label = F, group.by = "orig.ident", 
         do.return = TRUE, no.legend = F, #colors.use = singler.colors,
         pt.size = 1,label.size = 8 )+
    ggtitle("Original")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

save(Lung, file = "./data/Lung_20181026.Rda")
Iname = load("./data/Lung_20181026.Rda")

#======1.4 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../seurat_resources/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- HumanGenes(Lung,cc.genes[1:43])
g2m.genes <- HumanGenes(Lung,cc.genes[44:97])
# Assign Cell-Cycle Scores
Lung <- CellCycleScoring(object = Lung, s.genes = s.genes, g2m.genes = g2m.genes, 
                         set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = Lung, features.plot = HumanGenes(Lung,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
# regressing out the difference between the G2M and S phase scores
Lung@meta.data$CC.Difference <- Lung@meta.data$S.Score - Lung@meta.data$G2M.Score
# view cell cycle scores and phase assignments
head(x = Lung@meta.data)

#======1.5 Add project id =========================
orig.ident = Lung@meta.data$orig.ident
batch.effect = as.numeric(factor(orig.ident,levels = samples))
names(batch.effect) = rownames(Lung@meta.data)
Lung <- AddMetaData(object = Lung, metadata = batch.effect, col.name = "batch.effect")
table(Lung@meta.data$batch.effect)
head(x = Lung@meta.data)
#---------
projects = rep(NA,length(orig.ident))
projects[orig.ident %in% c("170511","CU7","CU-11-Proximal","CU-11-Distal",
                          "UNC-38","UNC-42")] = "ddSeq"
projects[orig.ident %in% c("UNC-44-Proximal","UNC-44-Distal",
                           "Ad-UNC-44-Terminal")] = "10X"
names(projects) = rownames(Lung@meta.data)
Lung <- AddMetaData(object = Lung, metadata = projects, col.name = "projects")
table(Lung@meta.data$projects)
head(x = Lung@meta.data)
#======1.6 vars.to.regress ScaleData =========================
features_threshold <- data.frame(c("nUMI","nGene","batch.effect","percent.mito","CC.Difference"),
                        c(10000,2000,2.0,0.05,0.05))
for(i in 1:nrow(features_threshold)){
    jpeg(paste0(path,"S2_",features_threshold[i,1],".jpeg"), units="in", width=10, height=7,res=600)
    P <- SingleFeaturePlot.1(Lung, feature = features_threshold[i,1],
                        threshold= features_threshold[i,2])
    print(P)
    dev.off()
}
Lung <- ScaleData(object = Lung, 
                  model.use = "linear", do.par=T, do.center = T, do.scale = T,
                  #vars.to.regress = c("nUMI","percent.mito","batch.effect","CC.Difference"),
                  display.progress = T)
#======1.7 Performing MNN-based correction =========================
#https://bioconductor.org/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html#4_performing_mnn-based_correction
set.seed(100)
original <- lapply(samples, function(x) Lung@scale.data[Lung@var.genes, 
                                                 (Lung@meta.data$orig.ident %in% x)])
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.order=T,
                                             approximate=TRUE)))
dim(mnn.out$corrected)
rownames(mnn.out$corrected) = Lung@cell.names
colnames(mnn.out$corrected) = paste0("MNN_",1:ncol(mnn.out$corrected))
#Storing a new MNN
Lung <- SetDimReduction(object = Lung, reduction.type = "MNN", slot = "cell.embeddings",
                       new.data = mnn.out$corrected)
Lung <- SetDimReduction(object = Lung, reduction.type = "MNN", slot = "key", 
                       new.data = "MNN_")
remove(original);GC()
Lung <- SetAllIdent(Lung,id = "orig.ident")

jpeg(paste0(path,"DimPlot.jpeg"), units="in", width=10, height=7,res=600)
DimPlot(object = Lung, reduction.use = "MNN", pt.size = 0.5)
dev.off()
#======1.7 unsupervised clustering based on MNN =========================
Lung <- RunPCA(object = Lung, pc.genes = Lung@var.genes, pcs.compute = 100, 
               do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = Lung)
PCElbowPlot(object = Lung, num.pc = 100)
PCHeatmap(Lung, pc.use = c(1:3, 25:30), cells.use = 500, do.balanced = TRUE)

DimElbowPlot.1(object = Lung, reduction.type = "MNN", 
             dims.plot = 50,slot = "cell.embeddings")


Lung <- RunTSNE(object = Lung, reduction.use = "MNN", dims.use = 1:50, 
                do.fast = TRUE, perplexity= 30)

Lung <- FindClusters(object = Lung, reduction.type = "MNN", 
                    dims.use = 1:50, resolution = 0.8, 
                     k.param = 30,force.recalc = T,
                     save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)

p2 <- TSNEPlot.1(object = Lung, do.label = F, group.by = "orig.ident", 
           do.return = TRUE, no.legend = T, 
           pt.size = 1,label.size = 4 )+
    ggtitle("Corrected")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"remove_batch.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1 +theme(legend.position="none"),p2)
dev.off()


p3 <- TSNEPlot.1(object = Lung, do.label = T, group.by = "ident", 
                 do.return = TRUE, no.legend = T, 
                 pt.size = 1,label.size = 6 )+
    ggtitle("Tsne plot for all clusters")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"tsneplot.jpeg"), units="in", width=10, height=7,res=600)
p3
dev.off()


p4 <- SplitTSNEPlot(object = Lung, split.by = "orig.ident",
                    do.label = F, group.by = "ident", 
                 do.return = TRUE, no.legend = T, 
                 pt.size = 1,label.size = 6 )

jpeg(paste0(path,"split_tsneplot.jpeg"), units="in", width=10, height=7,res=600)
print(do.call(plot_grid,p4))
dev.off()

save(Lung, file = "./data/Lung_20181101.Rda")
