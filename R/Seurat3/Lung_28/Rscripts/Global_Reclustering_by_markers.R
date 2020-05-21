########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","cowplot","magrittr","MAST",
                   "future","ggplot2","tidyr","harmony"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))

opts = data.frame(methods = rep(c("T20","IE", "DS", "3C","2607"),  time = 2),
                  labeled =  rep(c(T,  F), each = 5),
                  npcs = c(64,61,61,88,100,53,53,52,68,100),
                  stringsAsFactors = F)
(methods <- opts[i,])
(method = methods$methods)
(label = ifelse(methods$labeled, "labeled", "unlabeled"))
# ==================================================
step = 5

if(step == 1){
        # require 64GB
        # Find pc number 
        save.path <- paste0(path,method,"_",label,"/Find_pc_number/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        #======1.2 load  Seurat =========================
        df_samples <- readxl::read_excel("doc/Cell type markers for UMAP re-clustering.xlsx",
                                         sheet = sub("_.*","",methods$methods))
        object = readRDS(file = "data/Lung_28_Global_20200511.rds")
        DefaultAssay(object) = "SCT"
        Idents(object) = "annotations"
        if(!methods$labeled) object %<>% subset(idents = "unknown", invert = methods$labeled)
        Idents(object) = "orig.ident"
        object@neighbors = list()
        object@reductions = list()
        GC()
        # keep the VariableFeatures in marker list only
        marker_genes <- CaseMatch(search = df_samples$genes, match = rownames(object)) %>% as.character()
        
        VariableFeatures(object) = marker_genes
        print(length(VariableFeatures(object)))
        #====== find best pc number =========================
        object %<>% ScaleData
        object %<>% RunPCA(verbose = T,npcs = 100, features = VariableFeatures(object))

        npcs= 100
        
        p <- ElbowPlot(object, ndims = npcs)+
                ggtitle(paste("ElbowPlot for",methods))+
                TitleCenter()
        jpeg(paste0(save.path,"ElbowPlot.jpeg"),units="in", width=10, height=7,res=600)
        print(p)
        dev.off()
        
        a <- seq(1,97, by = 6)
        b <- a+5
        for(i in seq_along(a)){
                jpeg(paste0(save.path,"DimHeatmap_",i,"_",methods,"_",a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
                DimHeatmap(object, dims = a[i]:min(b[i],100),
                           nfeatures = 30,reduction = "pca")
                dev.off() 
        }
        
        object %<>% JackStraw(num.replicate = 20,dims = npcs)
        object %<>% ScoreJackStraw(dims = 1:npcs)
        a <- seq(1,100, by = 10)
        b <- a+9
        for(i in seq_along(a)){
                jpeg(paste0(save.path,"JackStrawPlot_",i,"_", a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
                print(JackStrawPlot(object, dims = a[i]:min(b[i],100)))
                Progress(i,length(a))
                dev.off()
        }
        saveRDS(object, file = paste0("output/Lung_28_",i,"-",method,"-",label,"_20200516.rds"))
}
# ReductionsPlots
if(step == 2){
        save.path <- paste0(path,method,"_",label,"/ReductionsPlots/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        
        object = readRDS(file = paste0("output/Lung_28_",i,"-",method,"-",label,"_20200516.rds"))
        
        ma <- function(x, n = 5) {stats::filter(x, rep(1 / n, n), sides = 1)}
        score.df <- JS(object = object[['pca']], slot = "overall")
        ma_Score = as.numeric(ma(score.df[,"Score"]))
        ma_Score[1:2] = ma_Score[3]
        print(paste("npcs =", npcs <- max(which(ma_Score < 0.05)+5))) # plus 5 just in case
        
        object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))
        
        jpeg(paste0(save.path,"RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
        system.time(object %<>% RunHarmony.1(group.by = "orig.ident", 
                                           dims.use = 1:npcs,
                                           tau = 0,
                                           theta = 2, plot_convergence = TRUE,
                                           nclust = 50, max.iter.cluster = 100))
        dev.off()
        
        object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
        object %<>% FindClusters(resolution = 0.8)
        object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
        object@assays$SCT@scale.data = matrix(0,0,0)
        format(object.size(object[["SCT"]]),unit = "GB")
        saveRDS(object, file = paste0("output/Lung_28_",i,"-",method,"-",label,"_20200516.rds"))
        
        Idents(object) = "orig.ident"
        lapply(c(TRUE, FALSE), function(lab) 
                UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.5,label = lab,
                    cols = Singler.colors,
                    label.repel = T,alpha = 0.9,
                    no.legend = T,label.size = 4, repel = T, title = "Harmony Integration",
                    do.print = T, do.return = F,save.path = save.path))
        
        Idents(object) = "cell_types"
        lapply(c(TRUE, FALSE), function(lab)
                UMAPPlot.1(object, group.by="cell_types",pt.size = 0.5,label = lab,
                    label.repel = T,alpha = 0.9,cols = ExtractMetaColor(object),
                    no.legend = T,label.size = 4, repel = T, title = "First annotation",
                    do.print = T, do.return = F,save.path = save.path))
        Idents(object) = "cell.labels"
        lapply(c(TRUE, FALSE), function(lab)
                UMAPPlot.1(object, group.by="cell.labels",pt.size = 0.5,label = lab,
                    label.repel = T,alpha = 0.9,cols = Singler.colors,
                    no.legend = T,label.size = 4, repel = T, title = "2nd annotation",
                    do.print = T, do.return = F,save.path = save.path))
        Idents(object) = "annotations"
        lapply(c(TRUE, FALSE), function(lab)
                UMAPPlot.1(object, group.by="annotations",pt.size = 0.5,label = lab,
                    cols = ExtractMetaColor(object),
                    label.repel = T, alpha= 0.9,
                    no.legend = T,label.size = 4, repel = T, title = "Last annotation",
                    do.print = T, do.return = F,save.path = save.path))
        FeaturePlot.1(object, features = c("MS4A1","CYTL1","TPSAB1",
                                           "DKK2","CHGA","L1CAM",
                                           "KRT14","MZB1","FOXI1"),
                      ncol = 3,do.print = T)
}

# serial resolution and generate seurat
if(step == 3){
        save.path <- paste0(path,method,"_",label,"/serial_resolutions/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        object = readRDS(file = paste0("output/Lung_28_",i,"-",method,"-",label,"_20200516.rds"))
        DefaultAssay(object) = "SCT"
        resolutions = c(seq(0.001,0.009, by = 0.001),seq(0.01,0.09, by = 0.01),seq(0.1,5, by = 0.1))
        for(i in 1:length(resolutions)){
                object %<>% FindClusters(resolution = resolutions[i])
                Idents(object) = paste0("SCT_snn_res.",resolutions[i])
                UMAPPlot.1(object, group.by=paste0("SCT_snn_res.",resolutions[i]),pt.size = 0.3,label = T,
                           label.repel = T,alpha = 0.9,
                           do.return = F,
                           no.legend = T,label.size = 4, repel = T, 
                           title = paste("res =",resolutions[i],"in",label,"_",method," based on harmony"),
                           do.print = T, save.path = save.path)
                Progress(i,length(resolutions))
        }
}

# Rshiny
if(step == 4){
        object = readRDS(file = paste0("output/Lung_28_",i,"-",method,"-",label,"_20200516.rds"))
        # by samples
        #(samples <- c("All_samples",sort(unique(object$orig.ident))))
        #Rshiny_path <- paste0("Rshiny/data/Lung_28_",method,"_",label,"_by_samples/")
        #PrepareShiny(object, samples, Rshiny_path, split.by = "orig.ident", 
        #             reduction = "umap",verbose = T)
        # by cell types
        (split.by = ifelse(i <= 4, "annotations","cell.labels") )
        object@meta.data[,split.by] %<>% as.character()
        (samples <- c("All_samples",sort(unique(as.character(object@meta.data[,split.by])))))
        Rshiny_path <- paste0("Rshiny/Lung_28_",method,"_",label,"_by_cell.types/")
        PrepareShiny(object, samples, Rshiny_path, split.by = split.by, 
                     reduction = "umap",verbose = T)
}
# re-run harmony
if(step == 5){
        # 64GB
        library(FrF2)
        set.seed(101)
        args = FrF2(16,4, factor.names=list(theta = c(0,8),
                                           sigma = c(0.1,1),
                                           nclust = c(20,100),
                                           tau = c(0,4)))
        (args= args[order(args$theta, args$sigma, args$nclust, args$tau),])
        rownames(args)=1:16
        (arg = args[i,])
        save.path <- paste0(path,"theta=",arg$theta,",sigma=",arg$sigma,",nclust=",arg$nclust,",tau=",arg$tau,"/Find_pc_number/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        #======1.2 load  Seurat =========================
        object = readRDS(file = "data/Lung_28_Global_20200511.rds")
        DefaultAssay(object) = "SCT"
        Idents(object) = "annotations"
        if(!methods$labeled) object %<>% subset(idents = "unknown", invert = methods$labeled)
        Idents(object) = "orig.ident"
        GC()
        # keep the VariableFeatures in marker list only
        df_samples <- readxl::read_excel("doc/Cell type markers for UMAP re-clustering.xlsx",
                                         sheet = sub("_.*","",methods$methods))
        marker_genes <- CaseMatch(search = df_samples$genes, match = rownames(object)) %>% as.character()
        
        VariableFeatures(object) = marker_genes
        print(length(VariableFeatures(object)))
        #===========
        object %<>% ScaleData
        (npcs = ncol(object@reductions$pca@cell.embeddings))
        print(paste("npcs =", npcs))
        object@reductions = list()
        object %<>% RunPCA(npcs = 100)
        jpeg(paste0(save.path,"RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
        system.time(object %<>% RunHarmony.1(group.by = "orig.ident", 
                                             dims.use = 1:npcs, sigma = arg$sigma,
                                             tau = arg$tau,
                                             theta = arg$theta, plot_convergence = TRUE,
                                             nclust = arg$nclust, max.iter.cluster = 100))
        dev.off()
        
        object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
        object %<>% FindClusters(resolution = 0.8)
        object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
        FeaturePlot.1(object, features = c("MS4A1","CYTL1","TPSAB1",
                                           "DKK2","CHGA","L1CAM",
                                           "KRT14","MZB1","FOXI1"),
                      ncol = 3,do.print = T, save.path = save.path)
}
# re-run harmony
if(step == 6){
        # 64GB
        #======1.2 load  Seurat =========================
        load(file = "data/Lung_28_20200116.Rda")
        object@assays$integrated =NULL
        object@assays$SCT@misc = NULL
        object@neighbors = list()
        
        Annotations = readRDS("Yang/proximal_distal_terminal_COPD/Subset_Reclustering_by_markers/Annotations.rds")
        Annotations = Unlist(Annotations)
        df = data.frame(annotations = names(Annotations), row.names = Annotations)
        unknown_cells <- colnames(object)[!(colnames(object) %in% Annotations)]
        df_unknown = data.frame(annotations = rep("unknown",length(unknown_cells)), row.names = unknown_cells)
        df %<>% rbind(df_unknown)
        df$barcode = rownames(df)
        df = df[colnames(object),]
        table(rownames(object@meta.data) == rownames(df))
        
        object[["annotations"]] = df$annotations
        object %<>% AddMetaColor(label= "annotations", colors = Singler.colors)
        Idents(object) = "annotations"
        if(!methods$labeled) object %<>% subset(idents = "unknown", invert = methods$labeled)
        Idents(object) = "orig.ident"
        GC()

        # keep the VariableFeatures in marker list only
        df_samples <- readxl::read_excel("doc/Cell type markers for UMAP re-clustering.xlsx",
                                         sheet = sub("_.*","",methods$methods))
        marker_genes <- CaseMatch(search = df_samples$genes, match = rownames(object)) %>% as.character()
        
        VariableFeatures(object) = marker_genes
        print(length(VariableFeatures(object)))
        #======= find npcs number ===============
        object %<>% ScaleData
        npcs <- 120
        object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))
        object %<>% JackStraw(num.replicate = 20,dims = npcs)
        object %<>% ScoreJackStraw(dims = 1:npcs)
        
        ma <- function(x, n = 5) {stats::filter(x, rep(1 / n, n), sides = 1)}
        score.df <- JS(object = object[['pca']], slot = "overall")
        ma_Score = as.numeric(ma(score.df[,"Score"]))
        ma_Score[1:2] = ma_Score[3]
        print(paste("npcs =", npcs <- max(which(ma_Score < 0.05)+5))) # plus 5 just in case
        
        object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))
        jpeg(paste0(save.path,"RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
        system.time(object %<>% RunHarmony.1(group.by = "orig.ident", 
                                             dims.use = 1:npcs,
                                             tau = 0,
                                             theta = 2, plot_convergence = TRUE,
                                             nclust = 100, max.iter.cluster = 100))
        dev.off()
        
        object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
        object %<>% FindClusters(resolution = 0.8)
        object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
        FeaturePlot.1(object, features = c("MS4A1","CYTL1","TPSAB1",
                                           "DKK2","CHGA","L1CAM",
                                           "KRT14","MZB1","FOXI1"),
                      ncol = 3,do.print = T, save.path = save.path)
        saveRDS(object, file = paste0("output/Lung_28_",i,"-",method,"-",label,"_20200516.rds"))
}
