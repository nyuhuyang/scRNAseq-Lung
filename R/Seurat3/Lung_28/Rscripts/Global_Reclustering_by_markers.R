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
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)


# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

opts = data.frame(methods = rep(c("T20","IE", "DS", "3C"),  time = 2),
                  labeled =  rep(c(T,  F), each = 4),
                  npcs = c(64,61,61,88,53,53,52,68),
                  stringsAsFactors = F)
(methods <- opts[args,])
(method = methods$methods)
(label = ifelse(methods$labeled, "labeled", "unlabeled"))
# ==================================================
step = 4

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
        object %<>% subset(idents = "unknown", invert = methods$labeled)
        Idents(object) = "orig.ident"
        object@neighbors = list()
        object@reductions = list()
        GC()
        # keep the VariableFeatures in marker list only
        marker_genes = FilterGenes(object, df_samples$genes)
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
        saveRDS(object, file = paste0("output/Lung_28_",args,"-",method,"-",label,"_20200516.rds"))
}
# ReductionsPlots
if(step == 2){
        save.path <- paste0(path,method,"_",label,"/ReductionsPlots/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        
        object = readRDS(file = paste0("output/Lung_28_",args,"-",method,"-",label,"_20200516.rds"))
        
        ma <- function(x, n = 5) {stats::filter(x, rep(1 / n, n), sides = 1)}
        score.df <- JS(object = object[['pca']], slot = "overall")
        ma_Score = as.numeric(ma(score.df[,"Score"]))
        ma_Score[1:2] = ma_Score[3]
        print(paste("npcs =", npcs <- max(which(ma_Score < 0.05)+5))) # plus 5 just in case
        
        object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))
        
        jpeg(paste0(save.path,"RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
        system.time(object %<>% RunHarmony(group.by.vars = "orig.ident", 
                                           assay.use="SCT",dims.use = 1:npcs,
                                           theta = 2, plot_convergence = TRUE,
                                           nclust = 50, max.iter.cluster = 100))
        dev.off()
        
        object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
        object %<>% FindClusters(resolution = 0.8)
        object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
        object@assays$SCT@scale.data = matrix(0,0,0)
        format(object.size(object[["SCT"]]),unit = "GB")
        saveRDS(object, file = paste0("output/Lung_28_",args,"-",method,"-",label,"_20200516.rds"))
        
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
}

# serial resolution and generate seurat
if(step == 3){
        save.path <- paste0(path,method,"_",label,"/serial_resolutions/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        object = readRDS(file = paste0("output/Lung_28_",args,"-",method,"-",label,"_20200516.rds"))
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
        object = readRDS(file = paste0("output/Lung_28_",args,"-",method,"-",label,"_20200516.rds"))
        # by samples
        #(samples <- c("All_samples",sort(unique(object$orig.ident))))
        #Rshiny_path <- paste0("Rshiny/data/Lung_28_",method,"_",label,"_by_samples/")
        #PrepareShiny(object, samples, Rshiny_path, split.by = "orig.ident", 
        #             reduction = "umap",verbose = T)
        # by cell types
        (split.by = ifelse(args <= 4, "annotations","cell.labels") )
        object@meta.data[,split.by] %<>% as.character()
        (samples <- c("All_samples",sort(unique(as.character(object@meta.data[,split.by])))))
        Rshiny_path <- paste0("Rshiny/Lung_28_",method,"_",label,"_by_cell.types/")
        PrepareShiny(object, samples, Rshiny_path, split.by = split.by, 
                     reduction = "umap",verbose = T)
}
# re-run harmony
if(step == 4){
        object = readRDS(file = paste0("output/Lung_28_",args,"-",method,"-",label,"_20200516.rds"))
        object %<>% ScaleData
        (npcs = ncol(object@reductions$pca@cell.embeddings))
        print(paste("npcs =", npcs))
        object[["harmony"]] = NULL
        object[["tsne"]] = NULL
        object[["umap"]] = NULL
        jpeg(paste0(save.path,"RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
        system.time(object %<>% RunHarmony.1(group.by = "orig.ident", 
                                           dims.use = 1:npcs,
                                           theta = 2, plot_convergence = TRUE,
                                           nclust = 50, max.iter.cluster = 100))
        dev.off()
        
        object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
        object %<>% FindClusters(resolution = 0.8)
        object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
        object@assays$SCT@scale.data = matrix(0,0,0)
        format(object.size(object[["SCT"]]),unit = "GB")
}