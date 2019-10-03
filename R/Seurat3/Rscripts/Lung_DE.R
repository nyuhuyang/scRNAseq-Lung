########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(gplots)
library(MAST)
library(ggpubr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
set.seed(1001)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# samples
samples = c("proximal","distal","terminal")
(con <- samples[args])
(load(file = paste0("data/Lung_24",con,"_20190918.Rda")))

if(con == "proximal"){
    print(unique(object@meta.data$conditions))
    Idents(object) = "integrated_snn_res.0.8"
    object %<>% RenameIdents("0" = "T cells:CD4+",
                             "1" = "T cells:TRM",
                             "2" = "Intermediate cells",
                             "3" = "Mucous gland cells",
                             "4" = "Goblet",
                             "5" = "Macrophages",
                             "6" = "B cells",
                             "7" = "T cells:TRM",
                             "8" = "Endothelial cells",
                             "9" = "Ciliated cells",
                             "10" = "Vascular Smooth muscle cells",
                             "11" = "T cells:7SK.2+",
                             "12" = "Endothelial cells",
                             "13" = "Basal cells",
                             "14" = "Fibroblasts",
                             "15" = "Neutrophil",
                             "16" = "Cartilage",
                             "17" = "hybrid cells",
                             "18" = "Mast cells",
                             "19" = "Myo-epithelial cells",
                             "20" = "Pericytes",
                             "21" = "NK cells",
                             "22" = "Proliferating basal cells",
                             "23" = "T cells:CD4+",
                             "24" = "Ionocytes/NEC",
                             "25" = "Neurons",
                             "26" = "Lymphatic endothelial cells",
                             "27" = "NK-like cells")
    object@meta.data$manual = as.character(Idents(object))
    object@meta.data = cbind(object@meta.data,object@reductions$umap@cell.embeddings)
    plasma <- object@meta.data$UMAP_1 >7 & object@meta.data$UMAP_2 < -10
    object@meta.data[plasma,"manual"] = "B cells:Plasma"
    pre_ciliated <- object$integrated_snn_res.0.8==9 & object@meta.data$UMAP_2 < 14.25 & object@meta.data$UMAP_1 > 1.5
    object@meta.data[pre_ciliated,"manual"] = "Pre-ciliated cells"
    basal_like <- object$integrated_snn_res.0.8==11 &  object@meta.data$UMAP_1 < 2.2 & object@meta.data$UMAP_2 > -2.2
    object@meta.data[basal_like,"manual"] = "Basal-like cells"
    ASM <- object$integrated_snn_res.0.8==19 &  object@meta.data$UMAP_1 < -7.8
    object@meta.data[ASM,"manual"] = "Airway smooth muscle"
    
    (load(file = "data/Epi_23proximal_20190904.Rda"))
    serous <- colnames(Epi)[Epi$integrated_snn_res.0.8 == 1]
    object@meta.data[serous,"manual"] = "Serous gland cells"
    ESML <- object$integrated_snn_res.0.8==10 &  object@meta.data$UMAP_2 < -5
    object@meta.data[ESML,"manual"] = "Endothelial smooth muscle-like cells"
    prlif_mf <- object$integrated_snn_res.0.8==5 &  object@meta.data$UMAP_1 < 7.2
    object@meta.data[prlif_mf,"manual"] = "Prolifereating macrophage"
    monocytes <- object$integrated_snn_res.0.8==5 &  object@meta.data$UMAP_2 < 2
    object@meta.data[monocytes,"manual"] = "Monocytes"
    DC <- object$integrated_snn_res.0.8==5 &  object@meta.data$UMAP_1 >10
    object@meta.data[DC,"manual"] = "Dendritic cells"
}

if(con == "distal"){
    DefaultAssay(object) <- 'integrated'
    for(i in c(15,16)/10){
        object %<>% FindClusters(resolution = i)
        Idents(object) = paste0("integrated_snn_res.",i)
    }
    DefaultAssay(object) = "RNA"
    print(unique(object@meta.data$conditions))
    Idents(object) = "integrated_snn_res.0.8"
    object %<>% RenameIdents("0" = "Ciliated cells",
                             "1" = "T cells:TRM",
                             "2" = "Fibroblasts",
                             "3" = "Endothelial cells:Capillary",
                             "4" = "Endothelial cells:HEVs",
                             "5" = "Alveolar type 2",
                             "6" = "Secretory cells",
                             "7" = "Smooth muscle",
                             "8" = "NK cells",
                             "9" = "Monocytes",
                             "10" = "Ciliated intermediate cells",
                             "11" = "Basal cells",
                             "12" = "Dendritic cells",
                             "13" = "B cells",
                             "14" = "Neutrophils",
                             "15" = "T cells:CD4+",
                             "16" = "Mast cells",
                             "17" = "Hybrid cells",
                             "18" = "Alveolar macrophages",
                             "19" = "Endothelial cells:Lymphatic",
                             "20" = "T cells:7SK.2+",
                             "21" = "Pre-ciliated cells",
                             "22" = "Endothelial cells",
                             "23" = "Plasma cells",
                             "24" = "Endothelial cells:Arterial",
                             "25" = "Neutrophils",
                             "26" = "Endothelial cells:Proliferating",
                             "27" = "Cartilage",
                             "28" = "Neuroendocrine")
                                 
    object@meta.data$manual = as.character(Idents(object))
    object@meta.data = cbind(object@meta.data,object@reductions$umap@cell.embeddings)
    object %<>% AddModuleScore(features = list(c("FCER1A","CCR7")), name = "DC")
    DC = object$DC1 > 3 & object$integrated_snn_res.0.8 %in% 12
    object@meta.data[DC,"manual"]="Dendritic cells"
    
    object@meta.data[object$integrated_snn_res.1.5 %in% 21,"manual"]="Secretory cells:Distal"
    object@meta.data[object$integrated_snn_res.1.5 %in% 28,"manual"]="Basal cells:Proliferating"
    Idents(object) = "integrated_snn_res.1.5"
    Proximal_markers <- read.csv("output/20190920/Proximal-2. DEGs for each proximal cell.types.csv",
                                 row.names = 1,stringsAsFactors = F) 
    Intermediate_genes <- Proximal_markers[(Proximal_markers$cluster %in% "Intermediate mucous-like"),"gene"]
    (Intermediate_genes %<>% head(5))
    Intermediate_genes <- FilterGenes(object,Intermediate_genes)
    object %<>% AddModuleScore(features = list(c(Intermediate_genes)), name = "Intermediate_genes")
    Intermediate = object$Intermediate_genes1 > 2 & object$integrated_snn_res.1.5 %in% 7 & object$UMAP_1 < -6
    object@meta.data[Intermediate,"manual"]="Intermediate cells"
    
    
    object@meta.data[object$integrated_snn_res.1.6 %in% 11,"manual"]="Smooth muscle:Vascular"
    object@meta.data[object$integrated_snn_res.1.6 %in% 21,"manual"]="Smooth muscle:Airway"
    object@meta.data[object$integrated_snn_res.1.6 %in% 27,"manual"]="Endothelial-smooth-muscle-like"
    object@meta.data[object$integrated_snn_res.1.6 %in% 28,"manual"]="Pericytes"
    object@meta.data[object$integrated_snn_res.1.6 %in% 31,"manual"]="Plasma cells"
    object@meta.data[object$integrated_snn_res.1.6 %in% 36,"manual"]="Alveolar type 1"
    object@meta.data[object$manual %in% "Endothelial cells","manual"] = "T cells:TRM"
    object@meta.data[object$manual %in% "Smooth muscle","manual"] = "Ciliated cells"
    }


if(con == "terminal"){
    print(unique(object@meta.data$conditions))
    object %<>% RenameIdents("0" = "T cells:CD8+",
                             "1" = "Alveolar macrophages",
                             "2" = "Alveolar type 2",
                             "3" = "Endothelial cells",
                             "4" = "Alveolar macrophages",
                             "5" = "NK cells",
                             "6" = "Alveolar macrophages",
                             "7" = "Ciliated cells",
                             "8" = "Distal secretory SCGB3A2+ SFTPB+",
                             "9" = "Smooth muscle",
                             "10" = "B cells", 
                             "11" = "Fibroblasts",
                             "12" = "Alveolar macrophages",
                             "13" = "Endothelial cells",
                             "14" = "Basal cells",
                             "15" = "T cells:CD8+",
                             "16" = "Alveolar macrophages",
                             "17" = "Mast cells",
                             "18" = "Endothelial cells",
                             "19" = "Endothelial cells",
                             "20" = "Lymphatic Endothelial cells",
                             "21" = "B cells:Plasma")
    object@meta.data$manual = as.character(Idents(object))
    object@meta.data = cbind(object@meta.data,object@reductions$umap@cell.embeddings)
    basal <- object@meta.data$manual %in% "Basal cells"
    hist(object@meta.data[basal,"UMAP_1"],breaks = 30)
    type1 <- object@meta.data$UMAP_1 < -1 & object@meta.data$UMAP_2 >6 
    object@meta.data[(basal & type1),"manual"] = "Alveolar type 1"
    
}
Idents(object) = "manual"
object %<>% sortIdent()
UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "manual", 
           cols = Singler.colors,legend.size = 15,
           do.return = F, no.legend = T, title = paste("UMAP plot for all clusters in",con),
           pt.size = 0.2,alpha = 0.85, label.size = 3, do.print = T,unique.name = "conditions")
UMAPPlot.1(object = object, label = F,label.repel = F, group.by = "manual", 
           cols = Singler.colors,legend.size = 15,
           do.return = F, no.legend = F, title = paste("UMAP plot for all clusters in",con),
           pt.size = 0.2,alpha = 0.85, label.size = 3, do.print = T,unique.name = "conditions")
UMAPPlot.1(object = object, label = F,label.repel = F, group.by = "manual", 
           cols = Singler.colors,legend.size = 15,
           do.return = F, no.legend = T, title = paste("UMAP plot for all clusters in",con),
           pt.size = 0.2,alpha = 0.85, label.size = 3, do.print = T,unique.name = "conditions")

Lung_markers <- FindAllMarkers.UMI(object, logfc.threshold = 0.1, only.pos = F,
                                   test.use = "MAST")
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
write.csv(Lung_markers,paste0(path,"Lung_24-",con,"_cell.types_markers.csv"))
save(object, file = paste0("data/Lung_24",con,"_20190918.Rda"))