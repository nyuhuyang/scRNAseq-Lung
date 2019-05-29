library(Seurat)
library(magrittr)
library(cowplot)
library(scran)
source("../R/Seurat_functions.R")
source("R/util.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#=== load data ======
(load(file="data/Lung_MNN_9_20181101.Rda"))
object <- SetAllIdent(object, id="orig.ident")
TSNEPlot.1(object,no.legend = F,do.label = T,do.print = T,do.return=F)

################################
# Bascal cells
################################
Bascal_cells <- FeaturePlot(object,"KRT15",do.identify = T)
write.csv(Bascal_cells,paste0(path,"Bascal_cells.csv"))

Bascal <- SubsetData(object, cells.use = Bascal_cells)
Bascal %<>% SetAllIdent(id="orig.ident")
TSNEPlot.1(Bascal,no.legend = F,do.label = T,do.print = T,do.return=F)

Bascal %<>% SetAllIdent(id="res.0.6")
p0 <- TSNEPlot.1(Bascal, do.print = T)
Bascal %<>% NormalizeData
Bascal %<>% FindVariableGenes(mean.function = ExpMean, 
                            dispersion.function = LogVMR, do.plot = F, 
                            x.low.cutoff = 0.125, x.high.cutoff = Inf, y.cutoff = 0.5)
length(Bascal@var.genes)
Bascal %<>% ScaleData %>% 
        RunPCA(pc.genes = Bascal@var.genes, pcs.compute = 30, do.print = F)
#====== Performing MNN-based correction =========================
#https://bioconductor.org/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html#4_performing_mnn-based_correction
set.seed(100)
(samples <- unique(Bascal@meta.data$orig.ident))
original <- lapply(samples, function(x) Bascal@scale.data[Bascal@var.genes, 
                                                          (Bascal@meta.data$orig.ident %in% x)])
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=30, auto.order=T,
                                             approximate=TRUE)))
dim(mnn.out$corrected)
rownames(mnn.out$corrected) = Bascal@cell.names
colnames(mnn.out$corrected) = paste0("MNN_",1:ncol(mnn.out$corrected))
#Storing a new MNN
Bascal <- SetDimReduction(object = Bascal, reduction.type = "MNN", slot = "cell.embeddings",
                                 new.data = mnn.out$corrected)
Bascal <- SetDimReduction(object = Bascal, reduction.type = "MNN", slot = "key", 
                                 new.data = "MNN_")
system.time({
        Bascal %<>% RunTSNE(reduction.use = "MNN", dims.use = 1:30, do.fast = TRUE) %>%
                FindClusters(reduction.type = "MNN", resolution = 0.6, dims.use = 1:30,
                             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                             force.recalc = TRUE, print.output = FALSE)
})

p3 <- TSNEPlot(Bascal, do.return = T, pt.size = 1, group.by = "orig.ident")
p4 <- TSNEPlot(Bascal, do.label = T, do.return = T, pt.size = 1)
jpeg(paste0(path,"S1_MMN_reTSNE.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Before re-run TSNE")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p4+ggtitle("After re-run TSNE")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
dev.off()


jpeg(paste0(path,"S1_MMN_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3+ggtitle("group by samples")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p4+ggtitle("group by clusters")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
dev.off()

FeaturePlot(Bascal,features.plot = "S100A2")
Bascal_markers <- FindAllMarkers.UMI(Bascal,logfc.threshold = 0.05,return.thresh = 0.05)
write.csv(Bascal_markers, paste0(path,"Bascal_markers.csv"))
h1 <- DoHeatmap.1(Bascal, Bascal_markers, Top_n = 50,
         title="Top 50 DE genes in each Bascal cells subpopulation",
         group.label.rot = F,cex.row = 4,remove.key =F,title.size = 12)
jpeg(paste0(path,"Bascal_cluster2.jpeg"), units="in", width=9.5, height=7,
     res=600)
print(h1)+ theme(strip.text.x = element_text(margin=margin(t = 30, r = 0, b = 0, l = 0)))
dev.off()

MakeHCorlorBar(Bascal, group_by = "res.0.6", remove.legend = F,
               file_name = "Bascal_cluster2",split.by = "orig.ident",do.print = TRUE,do.return=FALSE)
FeaturePlot(Bascal,features.plot =c("KRT5","S100A2","KRT14"),cols.use = c("grey", "red"))
g1 <- lapply(c("KRT5","S100A2","KRT14"), function(x) 
        SingleFeaturePlot.1(Bascal,feature = x, threshold = NULL))
jpeg(paste0(path,"Bascal1.jpeg"), units="in", width=10, height=3.5,res=600)
do.call(plot_grid,c(g1,nrow=1))
dev.off()

g2 <- lapply(c("MUC5B","SCGB1A1","MUC5AC"), function(x) 
        SingleFeaturePlot.1(Bascal,feature = x, threshold = NULL,gradient.use = c("grey", "green4")))
jpeg(paste0(path,"Bascal2.jpeg"), units="in", width=10, height=3.5,res=600)
do.call(plot_grid,c(g2,nrow=1))
dev.off()
################################
# Endothelial cells
################################
object %<>% SetAllIdent(id = "res.0.6")
TSNEPlot(object, do.label = T)
EC <- SubsetData(object,ident.use = 2)
keep <- FeaturePlot(EC, features.plot = "PECAM1",do.identify = T)
EC <- SubsetData(object,cells.use = keep)
p0 <- TSNEPlot.1(EC, do.print = T)
EC %<>% NormalizeData
EC %<>% FindVariableGenes(mean.function = ExpMean, 
                              dispersion.function = LogVMR, do.plot = F, 
                              x.low.cutoff = 0.125, x.high.cutoff = Inf, y.cutoff = 0.5)
length(EC@var.genes)
EC %<>% ScaleData %>% 
        RunPCA(pc.genes = EC@var.genes, pcs.compute = 30, do.print = F)
#====== Performing MNN-based correction =========================
#https://bioconductor.org/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html#4_performing_mnn-based_correction
set.seed(100)
(samples <- unique(EC@meta.data$orig.ident))
original <- lapply(samples, function(x) EC@scale.data[EC@var.genes, 
                                                          (EC@meta.data$orig.ident %in% x)])
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=30, auto.order=T,
                                             approximate=TRUE)))
dim(mnn.out$corrected)
rownames(mnn.out$corrected) = EC@cell.names
colnames(mnn.out$corrected) = paste0("MNN_",1:ncol(mnn.out$corrected))
#Storing a new MNN
EC <- SetDimReduction(object = EC, reduction.type = "MNN", slot = "cell.embeddings",
                          new.data = mnn.out$corrected)
EC <- SetDimReduction(object = EC, reduction.type = "MNN", slot = "key", 
                          new.data = "MNN_")
system.time({
        EC %<>% RunTSNE(reduction.use = "MNN", dims.use = 1:30, do.fast = TRUE) %>%
                FindClusters(reduction.type = "MNN", resolution = 0.4, dims.use = 1:30,
                             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                             force.recalc = TRUE, print.output = FALSE)
})

p3 <- TSNEPlot(EC, do.return = T, pt.size = 1, group.by = "orig.ident")
p4 <- TSNEPlot(EC, do.label = T, do.return = T, pt.size = 1)
jpeg(paste0(path,"S1_EC_reTSNE.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Before re-run TSNE")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p4+ggtitle("After re-run TSNE")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
dev.off()


jpeg(paste0(path,"S1_EC_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3+ggtitle("group by samples")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p4+ggtitle("group by clusters")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
dev.off()

EC_markers <- FindAllMarkers.UMI(EC,logfc.threshold = 0.05,return.thresh = 0.05)
write.csv(EC_markers, paste0(path,"EC_markers.csv"))
h2 <- DoHeatmap.1(EC, EC_markers, Top_n = 50,
                  title="Top 50 DE genes in each Endothelial cells subpopulation",
                  group.label.rot = F,cex.row = 4,remove.key =F,title.size = 12)
jpeg(paste0(path,"EC_cluster2.jpeg"), units="in", width=9.5, height=7,
     res=600)
print(h2)+ theme(strip.text.x = element_text(margin=margin(t = 30, r = 0, b = 0, l = 0)))
dev.off()

MakeHCorlorBar(EC, group_by = "res.0.4", remove.legend = F,
               file_name = "EC_cluster2",split.by = "orig.ident",do.print = TRUE,do.return=FALSE)
FeaturePlot(EC,features.plot =c("PECAM1","SELE"),cols.use = c("grey", "red"))

g3 <- lapply(c("PECAM1","SELE"), function(x) 
        SingleFeaturePlot.1(EC,feature = x, threshold = NULL))
jpeg(paste0(path,"EC_PECAM1_SELE.jpeg"), units="in", width=10, height=7,res=600)
do.call(plot_grid,c(g3,nrow=1))
dev.off()
