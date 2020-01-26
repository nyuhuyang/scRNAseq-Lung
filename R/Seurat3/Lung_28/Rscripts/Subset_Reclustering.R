########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   "magrittr","MAST","future"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# change the current plan to access parallelization
plan("multiprocess", workers = 4)
plan()

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

groups <- c("AT","B","D","En","F","Mon","SAE","SMG","SMP","T")
(g <- groups[args])

step = 5
# initial cluster
if(step == 1){
        (load(file = "data/Lung_28_20200116.Rda"))
        object$groups <- gsub("-.*","",object$cell_types) %>% gsub("[0-9]+","",.)
        object$groups %<>% plyr::mapvalues(from = c("AT","B","Car","D","DC",
                                                    "En","F","Mac","MC","Mon",
                                                    "NEC","Neu","Nr","P","PC",
                                                    "Per","SAE","SM","SMG","T"),
                                           to = c("AT","B","Car","D","Mon",
                                                  "En","F","Mon","MC","Mon",
                                                  "SAE","Neu","Nr","P","PC",
                                                  "SMP","SAE","SMP","SMG","T"))
        Idents(object) = "groups"
        object %<>% subset(idents = g)
        Idents(object) = "Doublets"
        object %<>% subset(idents = "Singlet")
        #======1.6 without intergration =========================
        DefaultAssay(object) = "SCT"
        object <- FindVariableFeatures(object = object, selection.method = "vst",
                                       num.bin = 20,nfeatures = 3000,
                                       mean.cutoff = c(0.1, 8), 
                                       dispersion.cutoff = c(1, Inf))
        object <- ScaleData(object = object,features = VariableFeatures(object))
        object <- RunICA(object, verbose =F,nics = 50)
        npcs =50
        object %<>% FindNeighbors(reduction = "ica",dims = 1:npcs)
        object %<>% FindClusters(resolution = 0.6)
        object %<>% RunTSNE(reduction = "ica", dims = 1:npcs, check_duplicates = FALSE)
        object %<>% RunUMAP(reduction = "ica", dims = 1:npcs)
        
        # cluster
        TSNEPlot.1(object, group.by = "SCT_snn_res.0.6",label = T,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   do.print = T,do.return = F, unique.name = "groups",
                   title = paste("Clusters in", g))
        UMAPPlot.1(object, group.by = "SCT_snn_res.0.6",label = T,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   do.print = T,do.return = F, unique.name = "groups",
                   title = paste("Clusters in", g))
        # cell types
        df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
        object$cell_types %<>% plyr::mapvalues(
                from = df_cell_types$`Cell types`,
                to = df_cell_types$Abbreviation)
        Idents(object) = "cell_types"
        object %<>% sortIdent()
        object <- AddMetaColor(object = object, label= "cell_types", colors = c(Singler.colors,Singler.colors))
        
        TSNEPlot.1(object, group.by = "cell_types",cols = ExtractMetaColor(object),label = T,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   do.print = T,do.return = F,
                   unique.name = "groups",
                   title = paste("Cell types in", g))
        UMAPPlot.1(object, group.by = "cell_types",cols = ExtractMetaColor(object),label = T,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   do.print = T,do.return = F,
                   unique.name = "groups",
                   title = paste("Cell types in", g))
        
        save(object, file = paste0("data/Lung_28_",g,"_20200121.Rda"))
}
# serial resolution and ICA heatmap
if(step == 2){
        path <- paste0(path,g,"/")
        if(!dir.exists(path))dir.create(path, recursive = T)
        load(file = paste0("data/Lung_28_",g,"_20200121.Rda"))
        DefaultAssay(object) <- 'SCT'
        
        res = c(1:12)/10
        for(i in seq_along(res)){
                object %<>% FindClusters(resolution = res[i])
                Idents(object) = paste0("SCT_snn_res.",res[i])
                UMAPPlot.1(object, group.by=paste0("SCT_snn_res.",res[i]),pt.size = 0.3,label = T,
                           label.repel = T,alpha = 0.9,
                           do.return = F,
                           no.legend = T,label.size = 4, repel = T, 
                           title = paste("res =",res[i],"in",g," based on ICA"),
                           unique.name = "groups",
                           do.print = T, save.path = path)
                Progress(i,length(res))
        }
        
        object <- RunICA(object, verbose =F,nics = 100)
        a <- seq(1,97, by = 6)
        b <- a+5
        for(i in seq_along(a)){
                jpeg(paste0(path,"DimHeatmap_ica_",g,"_",a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
                DimHeatmap(object, dims = a[i]:min(b[i],100),
                           nfeatures = 30,reduction = "ica")
                dev.off() 
        }
}

# test Integration
if(step == 3){
        path <- paste0(path,g,"/")
        if(!dir.exists(path))dir.create(path, recursive = T)
        load(file = paste0("data/Lung_28_",g,"_20200121.Rda"))
        #======1.6 intergration =========================
        DefaultAssay(object) = "SCT"
        Seurat_list <- SplitObject(object, split.by = "conditions")
        (size <- sapply(Seurat_list, ncol))
        Seurat_list[size < 20] = NULL
        Seurat_list %<>% lapply(SCTransform)
        object.features <- SelectIntegrationFeatures(Seurat_list, nfeatures = 3000)
        options(future.globals.maxSize= object.size(Seurat_list)*3)
        Seurat_list %<>% PrepSCTIntegration(anchor.features = object.features, verbose = FALSE)
        size <- sapply(Seurat_list, ncol)
        (k.filter <- min(size - 1, 200))
        (dims <- min(size - 2, 50))
        anchors <- FindIntegrationAnchors(Seurat_list, normalization.method = "SCT", 
                                          k.filter = k.filter,
                                          dims = 1:dims,
                                          anchor.features = object.features)
        remove(Seurat_list);GC()
        object <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight=30)
        remove(anchors);GC()
        object <- RunICA(object, verbose =F,nics = 50)
        npcs = 50
        object %<>% FindNeighbors(reduction = "ica",dims = 1:npcs)
        object %<>% FindClusters(resolution = 0.6)
        object %<>% RunTSNE(reduction = "ica", dims = 1:npcs, check_duplicates = FALSE)
        object %<>% RunUMAP(reduction = "ica", dims = 1:npcs)
        ElbowPlot()
        # cluster
        TSNEPlot.1(object, group.by = "integrated_snn_res.0.6",
                   label = T,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   do.print = T,do.return = F, unique.name = "groups",
                   save.path = path,
                   title = paste("Clusters in", g))
        UMAPPlot.1(object, group.by = "integrated_snn_res.0.6",
                   label = T,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   do.print = T,do.return = F, unique.name = "groups",
                   save.path = path,
                   title = paste("Clusters in", g))
        # cell types
        df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
        object$cell_types %<>% plyr::mapvalues(
                from = df_cell_types$`Cell types`,
                to = df_cell_types$Abbreviation)
        Idents(object) = "cell_types"
        object %<>% sortIdent()
        object <- AddMetaColor(object = object, label= "cell_types", colors = c(Singler.colors,Singler.colors))
        
        TSNEPlot.1(object, group.by = "cell_types",cols = ExtractMetaColor(object),label = T,
                   label.repel = T, pt.size = 1,label.size = 4, repel = T,no.legend = T,
                   do.print = T,do.return = F,
                   unique.name = "groups", save.path = path,
                   title = paste("Cell types in", g))
        UMAPPlot.1(object, group.by = "cell_types",cols = ExtractMetaColor(object),label = T,
                   label.repel = T, pt.size = 1,label.size = 4, repel = T,no.legend = T,
                   do.print = T,do.return = F,
                   unique.name = "groups", save.path = path,
                   title = paste("Cell types in", g))
        
        save(object, file = paste0("data/Lung_28_",g,"_20200122.Rda"))
}

# serial resolution and ICA heatmap on integrated data
if(step == 4){
        path <- paste0(path,g,"/")
        if(!dir.exists(path))dir.create(path, recursive = T)
        load(file = paste0("data/Lung_28_",g,"_20200122.Rda"))
        DefaultAssay(object) <- 'integrated'
        
        res = seq(0.01,0.09, by = 0.01)
        for(i in seq_along(res)){
                object %<>% FindClusters(resolution = res[i])
                Idents(object) = paste0("integrated_snn_res.",res[i])
                UMAPPlot.1(object, group.by=paste0("integrated_snn_res.",res[i]),pt.size = 0.3,label = T,
                           label.repel = T,alpha = 0.9,
                           do.return = F,
                           no.legend = T,label.size = 4, repel = T, 
                           title = paste("res =",res[i],"in",g," based on ICA"),
                           unique.name = "groups",
                           do.print = T, save.path = path)
                Progress(i,length(res))
        }
        
        object <- RunICA(object, verbose =F,nics = 100)
        a <- seq(1,50, by = 6)
        b <- a+5
        for(i in seq_along(a)){
                jpeg(paste0(path,"DimHeatmap_ica_",g,"_",a[i],"_",min(b[i],50),".jpeg"), units="in", width=10, height=7,res=600)
                DimHeatmap(object, dims = a[i]:min(b[i],50),
                           nfeatures = 30,reduction = "ica")
                dev.off() 
        }
}

if(step == 5){
        path <- paste0(path,g,"/")
        if(!dir.exists(path))dir.create(path, recursive = T)
        load(file = paste0("data/Lung_28_",g,"_20200122.Rda"))
        DefaultAssay(object) <- 'RNA'
        markers <- c("CDH5","CAPS","FOXJ1","DCN","LUM","KRT5",
                     "SFTPC","SRGN","LAPTM5","SCGB1A1", "CD3G","CSPG4")
        markers = FilterGenes(object, marker.genes = markers)
        FeaturePlot.1(object, markers,border = T,do.print = T, do.return = F,
                      unique.name = "groups", save.path = path)
}
FeaturePlot(object, markers)
# En (CDH5) 
# SAE:C (CAPS, FOXJ1) 
# F (DCN, LUM)
# SAE (KRT5, SFTPC,SCGB1A1)
# SMG (KRT5, SCGB1A1)
# SAE:Basal cells (KRT5) 
# AT (SFTPC,SCGB1A1+)
# B, T, Mon (SRGN, LAPTM5)
# SAE:Sec (SCGB1A1-very high)