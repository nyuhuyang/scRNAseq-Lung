########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#conda activate r4.0
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   "magrittr","harmony"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

########################################################################
#
#  1 read Seurat
# 
# ######################################################################
#======1.0 load  Seurat =========================
Lung_30 = readRDS(file = "data/Lung_30_20200702.rds")
(load(file = "data/Lung_8_time_20190808.Rda"))
Idents(object) ="orig.ident"
object %<>% subset(idents = c("Day-56","Day-122"), invert = T)

# ######################################################################
#1.1 – in vitro samples: d0, d3, d7, d14, d21, d28 combined.===================
save.path = paste0(path,"Lung_time_6/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
#1.2 run UMAP =========================

object %<>% ScaleData
npcs <- 105
object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))

object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

lapply(list(list(label = T, no.legend = T),list(label = F, no.legend = F),list(label = F, no.legend = T)), function(l) UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.1,
label = l$label,
           label.repel = T,alpha = 0.9, cols = Singler.colors,
           no.legend = l$no.legend, label.size = 4, repel = T, title = "No Integration with 6 time samples with 105 PCA",
           do.print = T, do.return = F, save.path = save.path))
format(object.size(object), unit = "GB")

data_df <- object@reductions$umap@cell.embeddings %>% as.data.frame() %>% 
    tibble::rownames_to_column("barcodes") %>% 
    left_join(object@meta.data %>% tibble::rownames_to_column("barcodes"), 
              by = "barcodes")
write.csv(data_df, file = paste0(save.path,basename(save.path),"_coordinates.csv"))

saveRDS(object, file = "data/Lung_time_6_20201001.rds")
#####################################################################################################

###################################################################################################
#2 – P, D, T, COPD samples P+D+T+COPD + in vitro d0, d3, d7, d14, d21, d28 combined.
#2.1 merge ===================================
save.path = paste0(path,"Lung_30_time_6/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

reference.list <- list("Lung_30" = Lung_30, "time" = object)
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), reference.list)
object@meta.data[,"batch"] = object@meta.data[,"conditions"]
object@meta.data[,"batch"] %<>% gsub("Day.*","Day",.)
object@meta.data[,"batch"] %<>% gsub("distal|proximal|terminal|COPD","region",.)

rm(reference.list);rm(Lung_30);GC()
remove(sce_list,sce_list2,Seurat_list);GC()

DefaultAssay(object)  = "SCT"
object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
#2.3 Find optimal PCA =========================
object %<>% ScaleData
npcs <- 150
object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))

jpeg(paste0(save.path,"S1_ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = npcs)
dev.off()

object %<>% JackStraw(num.replicate = 20,dims = 150)
object %<>% ScoreJackStraw(dims = 1:150)
a <- seq(1,150, by = 10)
b <- a+9
for(i in seq_along(a)){
    jpeg(paste0(save.path,"JackStrawPlot_",i,"_",a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a[i]:min(b[i],150)))
    Progress(i,length(a))
    dev.off()
}
#===
npcs <- 105
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "batch", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)

lapply(list(list(label = T, no.legend = T),
            list(label = F, no.legend = F),
            list(label = F, no.legend = T)), 
       function(l) UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.1,label = l$label,
           label.repel = T,alpha = 0.9, cols = Singler.colors,
           no.legend = l$no.legend, label.size = 4, repel = T, title = "No Integration with 36 samples with 105 PCA",
           do.print = T, do.return = F, save.path = save.path))

format(object.size(object), unit = "GB")
saveRDS(object, file = "data/Lung_30_time_6_20201009.rds")
#############################################################################################################################################


#############################################################################################################################################
#3 - epithelial P, D, T, COPD (epithelial cells you used in Monocle) + in vitro d0, d3, d7, d14, d21, d28.
save.path = paste0(path,"Epi_30_time_6/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

Idents(Lung_30) = "annotations3"
Epi_cells = c("BC","BC-S","BC-p","IC1","IC2",
               "IC-S","S","S-d","C1","C2","C3",
               "NEC","Ion","p-C","H","AT1","AT2",
               "AT2-1","AT2-p")
Lung_30 %<>% subset(idents = Epi_cells)
#3.1 merge ===================================
reference.list <- list("Lung_30" = Lung_30, "time" = object)
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), reference.list)
object@meta.data[,"batch"] = object@meta.data[,"conditions"]
object@meta.data[,"batch"] %<>% gsub("Day.*","Day",.)
object@meta.data[,"batch"] %<>% gsub("distal|proximal|terminal|COPD","region",.)

rm(reference.list);rm(Lung_30);GC()

#3.2 FindVariableFeatures =========================
DefaultAssay(object)  = "SCT"
object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
object %<>% ScaleData

#3.3 Find optimal PCA =========================
npcs <- 105
object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))

jpeg(paste0(save.path,"S1_ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = npcs)
dev.off()

jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "batch", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)

lapply(list(list(label = T, no.legend = T),list(label = F, no.legend = F),list(label = F, no.legend = T)), function(l) UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.1,
label = l$label,
           label.repel = T,alpha = 0.9, cols = Singler.colors,
           no.legend = l$no.legend, label.size = 4, repel = T, title = "No Integration with 36 Epi samples with 105 PCA",
           do.print = T, do.return = F, save.path = save.path))

format(object.size(object), unit = "GB")

object@assays$RNA@scale.data = matrix(0,0,0)
object@assays$SCT@scale.data = matrix(0,0,0)
object@assays$integrated = NULL
object@assays$SCT@misc = NULL
saveRDS(object, file = "data/Epi_30_time_6_20201009.rds")


# serial resolution and generate seurat
DefaultAssay(object) = "SCT"
resolutions = c(seq(0.001,0.009, by = 0.001),seq(0.01,0.09, by = 0.01),seq(0.1,5, by = 0.1))
for(i in 1:length(resolutions)){
        object %<>% FindClusters(resolution = resolutions[i])
        Idents(object) = paste0("SCT_snn_res.",resolutions[i])
        UMAPPlot.1(object, group.by=paste0("SCT_snn_res.",resolutions[i]),pt.size = 0.3,label = T,
                   label.repel = T,alpha = 0.9,
                   do.return = F,
                   no.legend = T,label.size = 4, repel = T, 
                   title = paste("res =",resolutions[i]),
                   do.print = T, save.path = paste0(save.path,"Serial_Resolution/"))
        Progress(i,length(resolutions))
}
