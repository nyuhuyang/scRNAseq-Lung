########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
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
                  stringsAsFactors = F)
(methods <- opts[args,])
method = methods$methods
(label = ifelse(methods$labeled, "labeled", "unlabeled"))
# ==================================================
step = 1

if(step == 1){
        # require 32GB
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
        
        object %<>% JackStraw(num.replicate = 20,dims = 100)
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
        # ReductionsPlots
        save.path <- paste0(path,method,"_",label,"/ReductionsPlots/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        
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
        object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs, check_duplicates = FALSE)
        object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
        
        Idents(object) = "orig.ident"
        lapply(c(TRUE, FALSE), function(lab) 
                UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.5,label = lab,
                    cols = ExtractMetaColor(object),
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
                fun(object, group.by="conditions",pt.size = 0.5,label = lab,
                    cols = ExtractMetaColor(object),
                    label.repel = T, alpha= 0.9,
                    no.legend = F,label.size = 4, repel = T, title = "Last annotation",
                    do.print = T, do.return = F,save.path = save.path))
        
        object@assays$RNA@scale.data = matrix(0,0,0)
        object@assays$SCT@scale.data = matrix(0,0,0)
        format(object.size(object[["RNA"]]),unit = "GB")
        saveRDS(object, file = paste0("output/Lung_28_",args,"-",method,"-",label,"_20200516.rds"))
}