########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################

library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("./data/")) dir.create("data")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/181230_Single_cell_sample list.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c(paste0("test",3)))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
conditions <- df_samples$conditions[sample_n]
projects <- df_samples$project[sample_n]
tests <- df_samples$tests[sample_n] 

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_3_20190109.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.seurat) %>%
        lapply(NormalizeData) %>%
        #lapply(ScaleData) %>%
        lapply(FindVariableGenes, do.plot = FALSE)

for(i in 1:length(samples)){
        object_list[[i]]@meta.data$tests <- tests[i]
        object_list[[i]]@meta.data$conditions <- conditions[i]
        object_list[[i]]@meta.data$projects <- projects[i]
        object_list[[i]]@meta.data$tissues <- tissues[i]
        
}
# we will take the union of the top 1k variable genes in each dataset for alignment
genes.use <- object_list %>% 
        lapply(function(object) head(rownames(object@hvg.info), 2000)) %>%
        unlist %>% unique
length(genes.use)

#========1.3 merge ===================================
object <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), object_list)
object@var.genes = genes.use
remove(sce_list,object_list);GC()

object = SetAllIdent(object, id = "orig.ident")
#======1.4 mito, QC, filteration =========================
(mito.genes <- grep(pattern = "^MT-", x = rownames(x = object@data), value = TRUE))
percent.mito <- Matrix::colSums(object@raw.data[mito.genes, ])/Matrix::colSums(object@raw.data)
object <- AddMetaData(object = object, metadata = percent.mito, col.name = "percent.mito")

(load(file = paste0("output/20190109/g1_3_20190109.Rda")))

object <- FilterCells(object = object, subset.names = c("nGene","nUMI","percent.mito"),
                   low.thresholds = c(500,200, -Inf), 
                   high.thresholds = c(Inf,Inf, 0.1))

object@ident = factor(object@ident,levels = samples)
g2 <- lapply(c("nGene", "nUMI", "percent.mito"), function(features){
        VlnPlot(object = object, features.plot = features, nCol = 3, 
                point.size.use = 0.2,size.x.use = 10, group.by = "ident",
                x.lab.rot = T, do.return = T)
})
save(g2,file= paste0(path,"g2_3_20190109.Rda"))

jpeg(paste0(path,"/S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                        scale_y_log10(limits = c(100,10000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                        scale_y_log10(limits = c(100,10000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                        scale_y_log10(limits = c(2000,100000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                        scale_y_log10(limits = c(2000,100000))))
dev.off()
jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                        ylim(c(0,0.25)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                        ylim(c(0,0.25))))
dev.off()
#======1.5 FindVariableGenes=======================
object <- NormalizeData(object = object)
jpeg(paste0(path,"/S1_dispersion.jpeg"), units="in", width=10, height=7,res=600)
object <- FindVariableGenes(object = object, mean.function = ExpMean, 
                            dispersion.function = LogVMR, do.plot = T, 
                            x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.4)
dev.off()
length(object@var.genes)
table(object@var.genes %in% genes.use)

#======1.5 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../R/seurat_resources/regev_lab_cell_cycle_genes.txt")
s.genes <- HumanGenes(object,cc.genes[1:43])
g2m.genes <- HumanGenes(object,cc.genes[44:97])
object <- CellCycleScoring(object = object, s.genes = s.genes, g2m.genes = g2m.genes, 
                        set.ident = FALSE)
jpeg(paste0(path,"/S1_RidgePlot.jpeg"), units="in", width=10, height=7,res=600)
RidgePlot(object = object, features.plot = HumanGenes(object,c("CCND1","CDK4","CCND2","CCND3")), #"CDK6",,"RB1"
          nCol = 2)
dev.off()
object@meta.data$CC.Difference <- object@meta.data$S.Score - object@meta.data$G2M.Score
object@meta.data$S.Score = object@meta.data$S.Score - min(object@meta.data$S.Score)
object@meta.data$G2M.Score = object@meta.data$G2M.Score - min(object@meta.data$G2M.Score)
tail(x = object@meta.data)


#======1.6 PCA =========================
object %<>% ScaleData %>%
        RunPCA(pc.genes = object@var.genes, pcs.compute = 30, do.print = F)
#Split.object <- SplitSeurat(object)
#Split.object %<>% lapply(ScaleData) %>%
#        lapply(function(object) RunPCA(object, pcs.compute = 30, do.print = F))

jpeg(paste0(path,"/S1_PCElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
PCElbowPlot(object, num.pc = 100)
dev.off()

jpeg(paste0(path,"/S1_PCHeatmap.jpeg"), units="in", width=10, height=7,res=600)
PCHeatmap(object, pc.use = c(1:3, 28:33), cells.use = 500, do.balanced = TRUE)
dev.off()

GC()
npca = 30
#rerun PCA

system.time({
        object %<>% RunTSNE(reduction.use = "pca", dims.use = 1:npca, do.fast = TRUE) %>%
                FindClusters(reduction.type = "pca", resolution = 0.6, dims.use = 1:npca,
                             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                             force.recalc = TRUE, print.output = FALSE)
})

#Split.object %<>% lapply(function(object) RunTSNE(object, reduction.use = "pca", dims.use = 1:npca, do.fast = TRUE)) %>%
#        lapply(function(object) FindClusters(object, reduction.type = "pca", resolution = 0.6, dims.use = 1:npca,
#                     save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
#                     force.recalc = TRUE, print.output = FALSE))

p0 <- DimPlot(object = object, reduction.use = "tsne", pt.size = 0.3, group.by = "orig.ident", do.return = T)
#======1.6 RunHarmony=======================
jpeg(paste0(path,"/S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony("orig.ident", dims.use = 1:30,
                                theta = 2, plot_convergence = TRUE,
                                nclust = 50, max.iter.cluster = 100))
dev.off()

object@ident %<>% factor(levels = samples)
p1 <- DimPlot(object = object, reduction.use = "harmony", pt.size = 0.3, group.by = "orig.ident", do.return = T)
p2 <- VlnPlot(object = object, features.plot = "Harmony1", group.by = "orig.ident", do.return = TRUE,
              x.lab.rot = T)
jpeg(paste0(path,"/S1_Harmony_vplot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1,p2)
dev.off()

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
system.time({
        object %<>% RunTSNE(reduction.use = "harmony", dims.use = 1:npca, do.fast = TRUE) %>%
                FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = 1:npca,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE)
})

p3 <- TSNEPlot(object, do.return = T, pt.size = 0.3, group.by = "orig.ident")
p4 <- TSNEPlot(object, do.label = T, do.return = T, pt.size = 0.3)

jpeg(paste0(path,"/S1_pca_vs_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Raw data")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p3+ggtitle("After alignment")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) )
dev.off()

jpeg(paste0(path,"/S1_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3+ggtitle("group by samples")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p4+ggtitle("group by clusters")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
dev.off()

(samples <- unique(object.all@meta.data$orig.ident) %>% sort)
for(i in 1:3){
        Split.object[[i]] <- SetAllIdent(Split.object[[i]], id = "res.0.6")
        
        g_Harmony <- TSNEPlot.1(object = Split.object[[i]], do.label = F, group.by = "ident",
                        do.return = TRUE, no.legend = F, 
                        #colors.use = ExtractMetaColor(Split.object[[i]]),
                        pt.size = 1,label.size = 4 )+
        ggtitle(paste("Tsne plot of all clusters in",
                      unique(Split.object[[i]]@meta.data$conditions)))+
        theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 
        
        jpeg(paste0(path,"TSNEplot-",i,".jpeg"), units="in", width=10, height=7,res=600)
        print(g_Harmony)
        dev.off()
}

object <- SetAllIdent(object, id = "res.0.6")
g_Harmony <- TSNEPlot.1(object = object, do.label = F, group.by = "ident",
                        do.return = TRUE, no.legend = F, 
                        colors.use = ExtractMetaColor(object),
                        pt.size = 1,label.size = 4 )+
        ggtitle(paste("Tsne plot of all clusters for all samples"))+
        theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot.jpeg"), units="in", width=10, height=7,res=600)
print(g_Harmony)
dev.off()


saveRDS(object@scale.data, file = "data/Lung.scale.data_Harmony_3_20190109.rds")
object@scale.data =NULL
save(object, file = "data/Lung_Harmony_3_20190109.Rda")


