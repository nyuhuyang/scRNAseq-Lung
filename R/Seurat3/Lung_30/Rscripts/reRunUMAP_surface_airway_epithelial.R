########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#conda activate r4.0
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

########################################################################
#
#  1 re-run UMAP for surface airway epithelial: COPD, D, P
# 
# ######################################################################
# To show that based on region-identity genes, COPD shifts from distal closer to proximal pattern.
# try cluster (UMAP, PCA) surface airway epithelial: COPD, D, P
# Input genes: DEGs: option 1 - D vs P; option 2 - D+T vs P
# criteria 1: all adj p <0.05
# criteria 2: adj p <0.05 and log2FC 0.5 or higher and -0.5 or lower
object = readRDS(file = "data/Lung_30_20200710.rds") 
DefaultAssay(object) = "SCT"

object$cell_types %<>% gsub("d-S","TASC",.)
SAE <- c("BC1","BC2","BC-p","IC1","IC2","IC3","S","TASC","H","p-C","C1","C2","C3","Ion","NE")
Idents(object) = "cell_types"
object %<>% subset(idents = SAE)
Idents(object) = "conditions"
Versions = list("with_T" = c("distal","terminal","COPD","proximal"),
                "without_T" = c("distal","COPD","proximal"))
version = Versions[args]

if(names(version) == "without_T") object %<>% subset(idents = version$without_T)

read.path <- "Yang/Lung_30/DE_analysis/surface_airway_epithelial/"
csv_list <- list("option_2" = "02_distal+terminal_vs_proximal",
                 "option_1" = "01_distal_vs_proximal")
csv = csv_list[[args]]
deg <- read.csv(paste0(read.path,"Lung_30_",csv,".csv"))
deg$pct.d = abs(deg$pct.1 - deg$pct.2)
print("table(abs(deg$avg_logFC) > 0.1 & deg$pct.d > 0.1)")
print(table(abs(deg$avg_logFC) > 0.1 & deg$pct.d > 0.1))
print("table(abs(deg$avg_logFC) > 0.5)")
print(table(abs(deg$avg_logFC) > 0.5))

#deg1[," -log10_p_adj"] = -log10(deg1$p_val_adj + 10^-150)
#deg1[,"log10_avg_UMI.1"] = log10(deg1$avg_UMI.1 + 1)
#deg1$abs_avg_logFC = abs(deg1$avg_logFC)
#deg1$pct.d = abs(deg1$pct.1 - deg1$pct.2)

#ggscatter(deg1, x = "avg_logFC", y = " -log10_p_adj", color = "pct.d")+
#        gradient_color(c("white","yellow","orange","chocolate1","red","darkred"))+
#        geom_vline(xintercept=0.15)+
#        geom_vline(xintercept=-0.1)

colnames(object[["umap"]]@cell.embeddings) %<>% paste0("orig_",.)
object[["orig.umap"]] <- CreateDimReducObject(embeddings = object[["umap"]]@cell.embeddings,
                                             key = "origUMAP_", assay = DefaultAssay(object))
#Option 1: abs(deg$avg_logFC) > 0.1 & deg$pct.d > 0.1, ~2.500 genes
object@reductions$umap = NULL
DefaultAssay(object) = "SCT"
object %<>% ScaleData(features = unique(deg[abs(deg$avg_logFC) > 0.1 & deg$pct.d > 0.1, "gene"]))
object %<>% RunPCA(npcs = 50, verbose = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:50)
colnames(object[["umap"]]@cell.embeddings) %<>% paste0("option1_",.)
object[["option1.umap"]] <- CreateDimReducObject(embeddings = object[["umap"]]@cell.embeddings,
                                              key = "option1UMAP_", assay = DefaultAssay(object))

#Option 2: abs(deg$avg_logFC) > 0.1 & deg$pct.d > 0.1, ~2.500 genes
object@reductions$umap = NULL
DefaultAssay(object) = "SCT"
object %<>% ScaleData(features = unique(deg[abs(deg$avg_logFC) > 0.5, "gene"]))
object %<>% RunPCA(npcs = 50, verbose = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:50)
colnames(object[["umap"]]@cell.embeddings) %<>% paste0("option2_",.)
object[["option2.umap"]] <- CreateDimReducObject(embeddings = object[["umap"]]@cell.embeddings,
                                                 key = "option2UMAP_", assay = DefaultAssay(object))
saveRDS(object, file = paste0("data/reRunMAP_SAE_",names(version),".rds"))

######################################################
# aggregate by sample
######################################################
object = readRDS(file = "data/Lung_30_20200710.rds") 

DefaultAssay(object) = "SCT"

object$cell_types %<>% gsub("d-S","TASC",.)
SAE <- c("BC1","BC2","BC-p","IC1","IC2","IC3","S","TASC","H","p-C","C1","C2","C3","Ion","NE")
Idents(object) = "cell_types"
object %<>% subset(idents = SAE)
Idents(object) = "orig.ident"
exp <- AverageExpression(object, assays = "SCT",slot = "data")

head(exp$SCT)

bulk <- CreateSeuratObject(counts = exp$SCT)
bulk$orig.ident = colnames(bulk)
meta.data = object@meta.data
meta.data = meta.data[!duplicated(meta.data$orig.ident),]
rownames(meta.data) = meta.data$orig.ident
meta.data = meta.data[rownames(bulk@meta.data),]
bulk@meta.data = meta.data

bulk %<>% FindVariableFeatures()
VariableFeatures(bulk) = unique(deg[abs(deg$avg_logFC) > 0.1 & deg$pct.d > 0.1, "gene"])
bulk %<>% ScaleData(features = VariableFeatures(bulk))
bulk %<>% RunPCA(npcs = 14, verbose = FALSE)
ElbowPlot(bulk)
Idents(bulk) = "conditions"
PCAPlot(bulk)
bulk %<>% RunUMAP(reduction = "pca", dims = 1:14)
colnames(bulk[["umap"]]@cell.embeddings) %<>% paste0("option1_",.)
bulk[["option1.umap"]] <- CreateDimReducObject(embeddings = bulk[["umap"]]@cell.embeddings,
                                                 key = "option1UMAP_", assay = DefaultAssay(bulk))
colnames(bulk[["pca"]]@cell.embeddings) %<>% paste0("option1_",.)
bulk[["option1.pca"]] <- CreateDimReducObject(embeddings = bulk[["pca"]]@cell.embeddings,
                                               key = "option1PCA_", assay = DefaultAssay(bulk))

#Option 2: abs(deg$avg_logFC) > 0.1 & deg$pct.d > 0.1, ~2.500 genes
bulk@reductions$umap = NULL
bulk@reductions$pca = NULL
VariableFeatures(bulk) = unique(deg[abs(deg$avg_logFC) > 0.5, "gene"])
bulk %<>% ScaleData(features = VariableFeatures(bulk))
bulk %<>% RunPCA(npcs = 14, verbose = FALSE)
ElbowPlot(bulk)
Idents(bulk) = "conditions"
PCAPlot(bulk)
bulk %<>% RunUMAP(reduction = "pca", dims = 1:14)
colnames(bulk[["umap"]]@cell.embeddings) %<>% paste0("option2_",.)
bulk[["option2.umap"]] <- CreateDimReducObject(embeddings = bulk[["umap"]]@cell.embeddings,
                                               key = "option2UMAP_", assay = DefaultAssay(bulk))
colnames(bulk[["pca"]]@cell.embeddings) %<>% paste0("option2_",.)
bulk[["option2.pca"]] <- CreateDimReducObject(embeddings = bulk[["pca"]]@cell.embeddings,
                                              key = "option2PCA_", assay = DefaultAssay(bulk))
bulk@reductions$umap = NULL
bulk@reductions$pca = NULL
saveRDS(bulk, file = paste0("data/Pseudobulk_SAE.rds"))

# volcano: COPD vs D; COPD vs P (comparing samples, not cells)
# Based on option 1 and option 2 genes separately. 
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")

bulk <- readRDS("data/Pseudobulk_SAE.rds")
Idents(bulk) = "conditions"
marker1 = FindMarkers_UMI(bulk,ident.1 = "COPD", ident.2 = "distal",logfc.threshold = 0,  min.pct = 0)
marker1$gene = rownames(marker1)
marker1 = marker1[!is.nan(marker1$p_val_adj),]
write.csv(marker1, paste0(path,"COPD_vs_D_SAE_bulk.csv"))
jpeg(paste0(path,"VolcanoPlots_COPD_vs_D_SAE_bulk.jpeg"),units="in", width=10, height=7,res=600)
VolcanoPlots(marker1,cut_off = "p_val",cut_off_logFC = 0.1,size=4) + ggtitle("COPD vs D in SAE Pseudobulk")
dev.off()

marker2 = FindMarkers_UMI(bulk,ident.1 = "COPD", ident.2 = "proximal",logfc.threshold = 0,  min.pct = 0)
marker2$gene = rownames(marker2)
marker2 = marker2[!is.nan(marker2$p_val_adj),]
write.csv(marker2, paste0(path,"COPD_vs_P_SAE_bulk.csv"))
jpeg(paste0(path,"VolcanoPlots_COPD_vs_P_SAE_bulk.jpeg"),units="in", width=10, height=7,res=600)
VolcanoPlots(marker2,cut_off = "p_val",cut_off_logFC = 0.1,size=4,top = 10) + ggtitle("COPD vs P in SAE Pseudobulk")
dev.off()

#average expression of these genes (option 1, option 2) per sample.
read.path <- "Yang/Lung_30/DE_analysis/surface_airway_epithelial/"
deg <- read.csv(paste0(read.path,"Lung_30_02_distal+terminal_vs_proximal.csv"))
deg$pct.d = abs(deg$pct.1 - deg$pct.2)

write.csv(bulk[["RNA"]]@data[deg[abs(deg$avg_logFC) > 0.1 & deg$pct.d > 0.1, "gene"],],
          file = paste0(path, "option1_avgExpr_by_sample.csv"))
write.csv(bulk[["RNA"]]@data[deg[abs(deg$avg_logFC) > 0.5, "gene"],],
          file = paste0(path, "option2_avgExpr_by_sample.csv"))

write.csv(exp$SCT, file = paste0(path, "SAE_expression.csv"))
