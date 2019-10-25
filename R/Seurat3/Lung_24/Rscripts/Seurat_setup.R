########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","magrittr",
                   "tidyr","gplots","MAST"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)


# SLURM_ARRAY_TASK_ID
set.seed=101
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# samples
samples = c("proximal","distal","terminal")
(con <- samples[args])

# Load Seurat
(load(file="data/Lung_24_20190918.Rda"))
Idents(object) = "conditions"
object %<>% subset(idents = con)
table(object$orig.ident)
#======1.5 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
npcs =100
object_list %<>% lapply(function(x) {
        x %<>% RunPCA(features = object.features, verbose = FALSE,npcs = npcs)
})
options(future.globals.maxSize= object.size(object_list)*3)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features, 
                                  verbose = FALSE)
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT", 
                                  anchor.features = object.features,
                                  reference = c(1, 2), reduction = "rpca", 
                                  dims = 1:npcs)
object <- IntegrateData(anchorset = anchors, normalization.method = "SCT",dims = 1:npcs)
save(object, file = paste0("data/Lung_24",con,"_20190918.Rda"))
remove(anchors,object_list);GC()
# FindNeighbors
object %<>% RunPCA(npcs = 100, verbose = FALSE)
npcs =100
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
UMAPPlot.1(object = object, group.by="integrated_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T,do.print = T,unique.name = "conditions",
           label.size = 4, repel = T,no.legend = T,
           title = paste("UMAP plot of",con, "samples"))

TSNEPlot.1(object, group.by="integrated_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T,do.print = T,unique.name = "conditions",
           label.size = 4, repel = T,no.legend = T,
           title = paste("tSNE plot of",con, "samples"))
object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = paste0("data/Lung_24",con,"_20190918.Rda"))


# Differential analysis
Idents(object) = "integrated_snn_res.0.8"
Lung_markers <- FindAllMarkers.UMI(object, logfc.threshold = 0.1, only.pos = T,
                                   test.use = "MAST")
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
write.csv(Lung_markers,paste0(path,"Lung_23-",con,"_markers.csv"))
Top_n = 5
top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
object %<>% ScaleData(features=unique(c(as.character(top$gene))))

DoHeatmap.1(object, marker_df = Lung_markers, Top_n = 5, do.print=T,
            angle = 0,group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
            assay = "SCT",
            label=T, cex.row=5, legend.size = 5,width=10, height=7,unique.name = T,
            title = paste("Top 5 markers of each clusters in",con,"sampels"))