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

# ========== cell types and markers ==============
# En (CDH5) 
# SAE:C (CAPS, FOXJ1) 
# F (DCN, LUM)
# SAE (KRT5, SFTPC,SCGB1A1)
# SMG (KRT5, SCGB1A1)
# SAE:Basal cells (KRT5) 
# AT (SFTPC,SCGB1A1+)
# B, T, Mon (SRGN, LAPTM5)
# SAE:Sec (SCGB1A1-very high)
step = 4
if(step == 1){
        # create folder
        save.path <- paste0(path,g,"/before-filter/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        # markers
        markers <- list(c("En","CDH5"),
                c("En","SRGN"),
                c("SAE:C","CAPS"),
                c("SAE:C","FOXJ1"),
                c("F","DCN"),
                c("F","LUM"),
                c("SAE","KRT5"),
                c("SAE","SFTPC"),
                c("SAE","SCGB1A1"),
                c("SMG","KRT5"),
                c("SMG","SCGB1A1"),
                #c("SMP","TAGLN"),
                c("SAE:BC","KRT5"),
                c("AT","SFTPC"),
                c("AT","SCGB1A1"),
                c("B","SRGN"),
                c("B","LAPTM5"),
                c("B","CD19"),
                c("T","SRGN"),
                c("T","LAPTM5"),
                c("T","CD3G"),
                c("Mon","SRGN"),
                c("Mon","LAPTM5"),
                c("SAE:Sec","SCGB1A1"))
        df_markers <- unlist(markers) %>% 
                matrix(nrow=length(markers), byrow=T) %>% 
                as.data.frame()
        colnames(df_markers) = c("cell_types", "markers")
        df_markers$cell_types %<>% gsub(":.*","",.)
        # load data 
        load(file = paste0("data/Lung_28_",g,"_20200122.Rda"))
        DefaultAssay(object) = "SCT"
        object@meta.data %<>% cbind(object@reductions$umap@cell.embeddings)
        # test features
        g_marker_genes <- df_markers[df_markers$cell_types %in% g,"markers"]
        all_marker_genes <- FilterGenes(object,unique(df_markers$markers))
        rm_genes = all_marker_genes[!all_marker_genes %in% g_marker_genes]
        FeaturePlot.1(object, features = all_marker_genes,
                      border = T, do.print = T, do.return = F,
                      unique.name = "groups",  save.path = save.path)
        object %<>% AddModuleScore(features = list(rm_genes),
                                   name = paste0(g,"_impurity"))
        fig <- FeaturePlot.1(object, features = paste0(g,"_impurity1"), 
                             title = paste(g,"impurity"), do.return = T)
        if(g == "AT") {
                fig = fig +
                geom_vline(xintercept = -6)+
                geom_vline(xintercept = 5)+
                geom_hline(yintercept = -7)+
                geom_hline(yintercept = 7)
                object %<>% subset(UMAP_1 > -6 & UMAP_1 < 5 &
                                   UMAP_2 > -7 & UMAP_2 < 7)
        }
        if(g == "B") {
                fig = fig +
                        geom_vline(xintercept = 1.5)+
                        geom_vline(xintercept = 6)+
                        geom_hline(yintercept = 4.5)
                object %<>% subset(UMAP_1 > 1.5 & UMAP_2 > 4.5, invert = T)
                object %<>% subset(UMAP_2 < 6)
        }
        if(g == "En") {
                fig = fig +
                        geom_vline(xintercept = -2)+
                        geom_vline(xintercept = 2.7)+
                        geom_vline(xintercept = 5)+
                        geom_hline(yintercept = 2.6)+
                        geom_hline(yintercept = 4.7)
                object %<>% subset(UMAP_1 > -2 & UMAP_2 > 4.7, invert = T)
                object %<>% subset(UMAP_1 > 2.7 & UMAP_2 > 2.6, invert = T)
                object %<>% subset(UMAP_1 < 5)
        }
        if(g == "F") {
                fig = fig +
                        geom_vline(xintercept = 5)+
                        geom_hline(yintercept = 6.5)+
                        geom_hline(yintercept = 1)
                object %<>% subset(UMAP_1 > 5 & UMAP_2 < 1, invert = T)
                object %<>% subset(UMAP_1 < 6.5)
        }
        if(g == "F") {
                fig = fig +
                        geom_vline(xintercept = 5)+
                        geom_hline(yintercept = 6.5)+
                        geom_hline(yintercept = 1)
                object %<>% subset(UMAP_1 > 5 & UMAP_2 < 1, invert = T)
                object %<>% subset(UMAP_1 < 6.5)
        }
        if(g == "Mon") {
                fig = fig +
                        geom_vline(xintercept = -9)+
                        geom_vline(xintercept = -2)+
                        geom_vline(xintercept = 4)+
                        geom_hline(yintercept = 6)+
                        geom_hline(yintercept = 1)
                object %<>% subset(UMAP_1 > -2 & UMAP_2 > 6, invert = T)
                object %<>% subset(UMAP_1 > 4 & UMAP_2 > 1, invert = T)
                object %<>% subset(UMAP_1 > -9)
        }
        if(g == "SAE") {
                fig = fig +
                        geom_vline(xintercept = -5)+
                        geom_vline(xintercept = 0)+
                        geom_hline(yintercept = -1)+
                        geom_hline(yintercept = 0)+
                        geom_hline(yintercept = 2.2)
                object %<>% subset(UMAP_1 < -5 & UMAP_2 < -1, invert = T)
                object %<>% subset(UMAP_1 > 0 & UMAP_2 > 0 & UMAP_2 < 2.2,
                                   invert = T)
        }
        if(g == "SMG") {
                fig = fig +
                        geom_vline(xintercept = -0.8)+
                        geom_vline(xintercept = 0.5)+
                        geom_hline(yintercept = 2.2)
                object %<>% subset(UMAP_1 > -0.8 & UMAP_1 < 0.5 &UMAP_2 > 2.2, invert = T)
        }
        if(g == "SMP") {
                fig = fig +
                        geom_vline(xintercept = -1.5)+
                        geom_vline(xintercept = 0)+
                        geom_hline(yintercept = -2.5)+
                        geom_hline(yintercept = -6)+
                        geom_hline(yintercept = -7)
                object %<>% subset(UMAP_1 > -1.5 & UMAP_2 > -6 & UMAP_2 < -2.5, invert = T)
                object %<>% subset(UMAP_1 > 0 & UMAP_2 > -7 & UMAP_2 < -2.5, invert = T)
        }
        if(g == "T") {
                fig = fig +
                        geom_vline(xintercept = -7)+
                        geom_vline(xintercept = -5)+
                        geom_vline(xintercept = -4)+
                        geom_hline(yintercept = 3.2)+
                        geom_hline(yintercept = 1.5)+
                        geom_hline(yintercept = 1)+
                        geom_hline(yintercept = -1)
                object %<>% subset(UMAP_1 < -5 & UMAP_2 > 3.2, invert = T)
                object %<>% subset(UMAP_1 < -7 & UMAP_2 > 1.5, invert = T)
                object %<>% subset(UMAP_1 > -7 & UMAP_1 < -5 & UMAP_2 < -1, invert = T)
                object %<>% subset(UMAP_1 > -5 & UMAP_1 < -4 & UMAP_2 < 1, invert = T)
        }
        jpeg(paste0(save.path,"UMAPPlot_",g,".jpeg"), 
             units="in", width=10, height=7,res=600)
        print(fig)
        dev.off()
        FeaturePlot.1(object, features = paste0(g,"_impurity1"),
                      title = paste(g,"impurity"), do.return = F,
                      do.print = T, unique.name = "groups",  
                      save.path = save.path)
        save(object, file = paste0("data/Lung_28_",g,"_20200126.Rda"))
}
if(step == 2){
        # create folder
        save.path <- paste0(path,g,"/re-run-umap/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        # load data
        load(file = paste0("data/Lung_28_",g,"_20200126.Rda"))
        DefaultAssay(object) = "SCT"
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
        object %<>% FindClusters(resolution = 0.1)
        object %<>% RunTSNE(reduction = "ica", dims = 1:npcs, check_duplicates = FALSE)
        object %<>% RunUMAP(reduction = "ica", dims = 1:npcs)
        # cluster
        UMAPPlot.1(object, group.by = "integrated_snn_res.0.1",
                   label = T,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   do.print = T,do.return = F, unique.name = "groups",
                   save.path = save.path,
                   title = paste("Clusters in", g))
        # cell types
        Idents(object) = "cell_types"
        object %<>% sortIdent()
        object <- AddMetaColor(object = object, label= "cell_types", colors = c(Singler.colors,Singler.colors))
        
        UMAPPlot.1(object, group.by = "cell_types",cols = ExtractMetaColor(object),label = T,
                   label.repel = T, pt.size = 1,label.size = 4, repel = T,no.legend = T,
                   do.print = T,do.return = F,
                   unique.name = "groups", save.path = save.path,
                   title = paste("Cell types in", g))
        
        save(object, file = paste0("data/Lung_28_",g,"_20200126.Rda"))
}

# serial resolution and ICA heatmap on integrated data
if(step == 3){
        res_path <- paste0(path,g,"/serial_resolution/")
        if(!dir.exists(res_path))dir.create(res_path, recursive = T)
        ica_path <- paste0(path,g,"/ica_heatmaps/")
        if(!dir.exists(ica_path))dir.create(ica_path, recursive = T)
        load(file = paste0("data/Lung_28_",g,"_20200126.Rda"))
        DefaultAssay(object) <- 'integrated'
        
        res = c(seq(0.01,0.1, by = 0.01),seq(0.2,1.6, by = 0.1))
        for(i in seq_along(res)){
                object %<>% FindClusters(resolution = res[i])
                Idents(object) = paste0("integrated_snn_res.",res[i])
                UMAPPlot.1(object, group.by=paste0("integrated_snn_res.",res[i]),pt.size = 0.3,label = T,
                           label.repel = T,alpha = 0.9,
                           do.return = F,
                           no.legend = T,label.size = 4, repel = T, 
                           title = paste("res =",res[i],"in",g," based on ICA"),
                           unique.name = "groups",
                           do.print = F, save.path = res_path)
                ElbowPlot.2(object, graph.name = "integrated_snn",)
                Progress(i,length(res))
        }
        
        object <- RunICA(object, verbose =F,nics = 100)
        a <- seq(1,50, by = 6)
        b <- a+5
        for(i in seq_along(a)){
                jpeg(paste0(ica_path,"DimHeatmap_ica_",
                            g,"_",a[i],"_",min(b[i],50),".jpeg"),
                     units="in", width=10, height=7,res=600)
                DimHeatmap(object, dims = a[i]:min(b[i],50),
                           nfeatures = 30,reduction = "ica")
                dev.off() 
        }
}

# generate ElbowPlot
if(step == 4){
        ElbowPlot_path <- paste0(path,g,"/ElbowPlot/")
        if(!dir.exists(ElbowPlot_path))dir.create(ElbowPlot_path, recursive = T)
        load(file = paste0("data/Lung_28_",g,"_20200126.Rda"))
        DefaultAssay(object) <- 'integrated'
        object <- RunICA(object, verbose =F,nics = 100)
        ElbowPlot.1(object, object, ndims = 100, reduction = "ica",
                    title = paste("Standard Deviation against ICA numbers in cell type",g),
                    do.print = T,do.return = F, save.path = ElbowPlot_path)
        res = c(seq(0.01,0.1, by = 0.01),seq(0.2,1.6, by = 0.1))
        for(i in seq_along(res)){
                object %<>% FindClusters(resolution = res[i])
                Progress(i,length(res))
        }
        ElbowPlot.2(object, graph.name = "integrated_snn",check.by = "Cluster.numbers",
                    unique.name = "groups",
                    title = paste("Total distance against different cluster numbers in cell type",g),
                    do.print = T,do.return = F, save.path = ElbowPlot_path)
        ElbowPlot.2(object, graph.name = "integrated_snn",check.by = "Resolutions",
                    unique.name = "groups",
                    title = paste("Total distance against different resolutions in cell type",g),
                    do.print = T,do.return = F, save.path = ElbowPlot_path)
}