########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(magrittr)
library(DoubletFinder)
library(kableExtra)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# change the current plan to access parallelization
########################################################################
#
#  2. DoubletFinder 
# 
# ######################################################################
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
object@assays$SCT = NULL

Idents(object) = "orig.ident"
(samples = unique(Idents(object)))
object_list <- SplitObject(object,split.by = "orig.ident")
GC()
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
npcs =100
sweep.res_list <- list()
for (i in 1:length(object_list)) {
    sweep.res_list[[i]] <- paramSweep_v4(object_list[[i]], PCs = 1:npcs, sct = T)
    message(paste("Complete", i,":", length(object_list)," ==================="))
}
save(sweep.res_list,file = "output/",con,"_sweep.res_list.Rda")

sweep_list <- lapply(sweep.res_list, function(x) summarizeSweep(x, GT = FALSE))
bcmvn_list <- lapply(sweep_list,find.pK)
# find histgram local maximam
find.localMaxima <- function(x) {
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
        y <- y[-1]
    }
    which(x == max(x[y]))
}

(maximal_pk <- sapply(bcmvn_list,function(x) {
    as.numeric(as.character(x[find.localMaxima(x$BCmetric),"pK"]))
    }))
maximal_pk

# http://rstudio-pubs-static.s3.amazonaws.com/329613_f53e84d1a18840d5a1df55efb90739d9.html
qplot_2axis <- function(data,x = "pK", y1 = "MeanBC", y2 = "BCmetric"){
    if(class(data[,x]) == "factor") data[,x] <- as.numeric(as.character(data[,x]))
    data_y1 <- data[,y1]
    data_y2 <- data[,y2]
    a <- range(data_y1)
    b <- range(data_y2)
    scale_factor <- diff(a)/diff(b)
    data_y2 <- ((data_y2 - b[1]) * scale_factor) + a[1]
    trans <- ~ ((. - a[1]) / scale_factor) + b[1]
    
    g <- ggplot(data = data, aes_string(x = x, y = y1))+
        geom_line()+geom_point()+
        geom_point(aes(y = data_y2),colour = "blue")+
        geom_line(aes(y = data_y2),colour = "blue")+
        scale_y_continuous(name = y1,
                           sec.axis = sec_axis(trans=trans, name=y2))+
        theme(axis.text.y.right = element_text(color = "blue"))
    
    g
    
}
qplot_2axis(data = bcmvn_list[[2]])

Multiplet_Rate <- function(object, numBatches = 1, num10xRuns = 1){
    
    numCellsRecovered = 1.0 * ncol(object)
    m = 4.597701e-06
    r = 0.5714286
    
    numCellsLoaded = numCellsRecovered / r
    multipletRate = m * numCellsLoaded / num10xRuns
    
    singletRate = 1.0 - multipletRate;
    numSinglet = singletRate * numCellsRecovered
    numMultiplet = numCellsRecovered - numSinglet
    numIdentMultiplet = numMultiplet * (numBatches - 1) / numBatches
    numNonIdentMultiplet = numMultiplet - numIdentMultiplet
    numCells = numSinglet + numNonIdentMultiplet
    
    return(numNonIdentMultiplet/numCells)
}
Multiplet_Rate(object_list[[1]])
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
for(i in 1:length(object_list)){
    print(paste("processing",unique(object_list[[i]]$orig.ident)))
    homotypic.prop <- modelHomotypic(object_list[[i]]@meta.data$cell.types)
    nExp_poi <- round(Multiplet_Rate(object_list[[i]])*length(colnames(object_list[[i]])))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    object_list[[i]] <- doubletFinder_v3(object_list[[i]], PCs = 1:50, 
                                                     pN = 0.25, pK = maximal_pk[i], 
                                                     nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    object_list[[i]] <- doubletFinder_v3(object_list[[i]], PCs = 1:50, pN = 0.25, maximal_pk[i],
                                   nExp = nExp_poi.adj, reuse.pANN = grep("pANN",colnames(object_list[[i]]@meta.data),value = T),
                                   sct = TRUE)
    colName = colnames(object_list[[i]]@meta.data)
    colName[grep("DF.classifications",colName)] = c("Low_confident_doublets",
                                                    "High_confident_doublets")
    colnames(object_list[[i]]@meta.data) = colName
}
for(i in 1:length(object_list)){
    object_list[[i]]@meta.data$row.names = rownames(object_list[[i]]@meta.data)
}
meta.data_list <- lapply(object_list, function(x) x@meta.data)
meta.data = bind_rows(meta.data_list)
rownames(meta.data) = meta.data$row.names
meta.data = meta.data[colnames(object),]
meta.data$doublets = gsub("Doublet","Doublet-Low Confidence",meta.data$Low_confident_doublets)
meta.data[meta.data$High_confident_doublets %in% "Doublet","doublets"] = "Doublet-High Confidence "
meta.data = cbind(object@meta.data,meta.data$doublets)
colnames(meta.data)[ncol(meta.data)] = "Doublets"

(load(file = paste0("data/Lung_24",con,"_20190918.Rda")))
object@meta.data = meta.data

TSNEPlot.1(object, group.by = "Doublets",cols = c("red","black"), 
           title = "Singlets and possible Doublets", do.print = T,pt.size = 0.3)
UMAPPlot.1(object, group.by = "Doublets",cols = c("red","black"), 
           title = "Singlets and possible Doublets", do.print = T,pt.size = 0.3)

table(object$Doublets) %>% prop.table %>% kable %>% kable_styling()
table(object$Doublets, object$cell.types) %>% kable %>% kable_styling()

save(object,file=paste0("data/Lung_24",con,"_20191002.Rda"))
