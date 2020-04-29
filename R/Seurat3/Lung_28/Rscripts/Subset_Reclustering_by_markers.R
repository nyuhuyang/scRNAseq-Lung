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

cell_types <- rep(c("T","En","C","Mon","BC_S",
                "F","AT","SM","SMG","Cr_Nr",
                "SAE","SAE_AT","SAE_AT_T","En_SM","F_Cr_Nr",
                "B_PC_Mon","Epithelial","Lymphoid","Mesenchymal","Myeloid"), each = 4)

methods = rep(c("T20","IE", "DS", "3C"), times = 20)
if(args > 80 & args <= 84){
        i = args; args=66
        resolution = c( 0.5,1.0,2.0,3.3)[i-80]
}
if(args > 84 & args <= 86){
        i = args; args=77
        resolution = c( 0.5,0.7)[i-84]
}
(method <- methods[args])
(cell_type = cell_types[args])
# ==================================================
step = 3
# Find pc number and generate seurat
if(step == 1){
        save.path <- paste0(path,args,"_",cell_type,"_",method,"/Find_pc_number/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        #======1.2 load  Seurat =========================
        df_samples <- readxl::read_excel("doc/Cell type markers for UMAP re-clustering.xlsx",
                                         sheet = sub("_.*","",method))
        object = readRDS(file = paste0("data/Lung_28_",ceiling(args/4),"-",cell_type,"_20200206.rds"))
        DefaultAssay(object) = "SCT"
        object@neighbors = list()
        object@reductions = list()
        GC()
        # keep the VariableFeatures in marker list only
        marker_genes = FilterGenes(object, df_samples$genes)
        VariableFeatures(object) = marker_genes
        #====== find best pc number =========================
        object %<>% ScaleData
        object %<>% RunPCA(verbose = T,npcs = 100, features = VariableFeatures(object))
        PCAPlot.1(object, do.print = T, save.path = save.path)
        
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
        saveRDS(object, file = paste0(path, "Lung_28_",args,"_",cell_type,"_",method,"-2020406.rds"))
        
        save.path <- paste0(path, args,"_",cell_type,"_",method,"/ReductionsPlots/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        #        object = readRDS(file = paste0("output/Lung_28_",args,"_",cell_type,"_",method,"-2020406.rds"))
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
        lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun) 
                fun(object, group.by="orig.ident",pt.size = 0.5,label = F,
                    cols = ExtractMetaColor(object),
                    label.repel = T,alpha = 0.9,
                    no.legend = T,label.size = 4, repel = T, title = "Harmony Integration",
                    do.print = T, do.return = F,save.path = save.path))
        
        meta.data = readRDS(file = "output/20200131/cell_types.rds")
        meta.data$barcode = rownames(meta.data)
        object1 = readRDS(file = "data/Lung_28_Global_20200219.rds") 
        object1$barcode = colnames(object1)
        meta.data %<>% full_join(object1@meta.data[,c("barcode","cell.labels",
                                                      "cell.labels.colors")],
                                 by = "barcode")
        rownames(meta.data) = meta.data$barcode
        meta.data = meta.data[rownames(object@meta.data),]
        object[["cell_types"]] = meta.data$cell_types
        object[["cell_types.colors"]] = meta.data$cell_types.colors  
        object[["cell.labels"]] = meta.data$cell.labels
        object[["cell.labels.colors"]] = meta.data$cell.labels.colors        
        Idents(object) = "cell_types"
        lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
                fun(object, group.by="cell_types",pt.size = 0.5,label = T,
                    label.repel = T,alpha = 0.9,cols = ExtractMetaColor(object),
                    no.legend = T,label.size = 4, repel = T, title = "First annotation",
                    do.print = T, do.return = F,save.path = save.path))
        Idents(object) = "cell.labels"
        lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
                fun(object, group.by="cell.labels",pt.size = 0.5,label = T,
                    label.repel = T,alpha = 0.9,cols = Singler.colors,
                    no.legend = T,label.size = 4, repel = T, title = "Final annotation",
                    do.print = T, do.return = F,save.path = save.path))
        
        lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
                fun(object, group.by="conditions",pt.size = 0.5,label = T,
                    cols = c('#601b3f','#3D9970','#FF4136','#FF851B'),
                    label.repel = T, alpha= 0.9,
                    no.legend = F,label.size = 4, repel = T, title = "Conditions",
                    do.print = T, do.return = F,save.path = save.path))
        
        object@assays$RNA@scale.data = matrix(0,0,0)
        object@assays$SCT@scale.data = matrix(0,0,0)
        format(object.size(object[["RNA"]]),unit = "GB")
        saveRDS(object, file = paste0("output/Lung_28_",args,"_",cell_type,"_",method,"-2020406.rds"))
}
# Find pc number and generate seurat
if(step == 2){
        object = readRDS(file = paste0("output/Lung_28_",args,"_",cell_type,"_",method,"-2020406.rds"))
        DefaultAssay(object) = "SCT"
        resolutions = c(seq(0.001,0.009, by = 0.001),seq(0.01,0.09, by = 0.01),seq(0.1,5, by = 0.1))
        for(i in 1:length(resolutions)){
                object %<>% FindClusters(resolution = resolutions[i])
                Idents(object) = paste0("SCT_snn_res.",resolutions[i])
                UMAPPlot.1(object, group.by=paste0("SCT_snn_res.",resolutions[i]),pt.size = 0.3,label = T,
                           label.repel = T,alpha = 0.9,
                           do.return = F,
                           no.legend = T,label.size = 4, repel = T, 
                           title = paste("res =",resolutions[i],"in",cell_type,"_",method," based on harmony"),
                           do.print = T, save.path = paste0(path,args,"_",cell_type,"_",method,"/"))
                Progress(i,length(resolutions))
        }
}

# DE analysis
if(step == 3){
        object = readRDS(file = paste0("output/Lung_28_",args,"_",cell_type,"_",method,"-2020406.rds"))
        DefaultAssay(object) = "SCT"
        object %<>% FindClusters(resolution = resolution)
        Idents(object) = paste0("SCT_snn_res.",resolution)
        res = FindAllMarkers.UMI(object,logfc.threshold = 0.5)
        write.csv(res, file = paste0(path,"DE_",args,"_",cell_type,"_",method,"_res=",resolution,".csv"))
}