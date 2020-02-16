########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
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

# re-run the KNN with and without intergration
reductions = c("harmony", "pca")
(reduction <- reductions[args])

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

ambigous = read.csv("Yang/proximal_distal_terminal_COPD/Harmony/Annotations/Lung_28_ambigous_nolabel_barcodes.csv",row.names = 1,
                    stringsAsFactors = F) %>% pull
g_name = "ambigous"
args = 21
# use ElbowPlot, PCA heatmap and JackStrawPlot to select PCA number 
if(step == 1){ # DimHeatmap 32GB, JackStraw need 128 GB
        save.path <- paste0(path,"21-ambigous","/Select_npcs/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        (load(file = "data/Lung_28_harmony_rmD_20200205.Rda"))
        object %<>% subset(cells = ambigous)
        GC()
        table(object$cell_types)
        DefaultAssay(object) = "SCT"
        object <- FindVariableFeatures(object = object, selection.method = "vst",
                                       num.bin = 20,nfeatures = 2000,
                                       mean.cutoff = c(0.1, 8), 
                                       dispersion.cutoff = c(1, Inf))
        object <- ScaleData(object = object,features = VariableFeatures(object))
        object <- RunPCA(object, verbose =F, npcs = 100)
        npcs = 100
        
        p <- ElbowPlot(object, ndims = npcs)+
                ggtitle(paste("ElbowPlot for",g_name))+
                TitleCenter()
        jpeg(paste0(save.path,"ElbowPlot.jpeg"),units="in", width=10, height=7,res=600)
        print(p)
        dev.off()
        
        a <- seq(1,97, by = 6)
        b <- a+5
        for(i in seq_along(a)){
                jpeg(paste0(save.path,"DimHeatmap_",i,"_",g_name,"_",a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
                DimHeatmap(object, dims = a[i]:min(b[i],100),
                           nfeatures = 30,reduction = "pca")
                dev.off()
                Process(i,length(a))
        }
        
        object %<>% JackStraw(num.replicate = 20,dims = 100)
        object %<>% ScoreJackStraw(dims = 1:npcs)
        a <- seq(1,100, by = 10)
        b <- a+9
        for(i in seq_along(a)){
                jpeg(paste0(save.path,"JackStrawPlot_",i,"_", g_name, a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
                print(JackStrawPlot(object, dims = a[i]:min(b[i],100)))
                Progress(i,length(a))
                dev.off()
        }
        saveRDS(object, file = paste0("data/Lung_28_21-",g_name,"_20200206.rds"))
}

# run harmony
if(step == 2){ # need 16 GB ?
        save.path <- paste0(path,"21-",g_name,"/subset_rerun_Harmony/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        object = readRDS(file = paste0("data/Lung_28_21-",g_name,"_20200206.rds"))
        #======1.6 intergration =========================
        (npcs = 90)
        object <- RunPCA(object, verbose =F, npcs = npcs)
        if(reduction = "harmony"){
                jpeg(paste0(save.path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
                system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                                     theta = 2, plot_convergence = TRUE,
                                                     nclust = 50, max.iter.cluster = 100))
                dev.off()
        }

        object %<>% FindNeighbors(reduction = reduction,dims = 1:npcs)
        object %<>% RunTSNE(reduction = reduction, dims = 1:npcs, check_duplicates = FALSE)
        object %<>% RunUMAP(reduction = reduction, dims = 1:npcs)
        saveRDS(object, file = paste0("data/Lung_28_21-",g_name,"-",reduction,"_20200206.rds"))
        

        UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.5,label = F,
            label.repel = T,alpha = 0.9,
            cols = Singler.colors,
            no.legend = T,label.size = 4, repel = T, 
            title = paste("Harmony Integration in",g_name),
            do.print = T, do.return = F,save.path = save.path)
        UMAPPlot.1(object, group.by="conditions",pt.size = 0.5,label = T,
            cols = c('#601b3f','#3D9970','#FF4136','#FF851B'),
            label.repel = T, alpha= 0.9,
            no.legend = F,label.size = 4, repel = T, 
            title = paste("Harmony Integration in",g_name),
            do.print = T, do.return = F,save.path = save.path)
}

if(step == 3){        # test serial_resolutions
        save.path <- paste0(path,args,"-",g_name,"/serial_resolutions/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        object = readRDS(file = paste0("data/Lung_28_",args,"-",g_name,"_20200206.rds"))
        resolutions = c(seq(0.001,0.009, by = 0.001),seq(0.01,0.09, by = 0.01),seq(0.1,3, by = 0.1))
        for(i in seq_along(resolutions)){
                object %<>% FindClusters(resolution = resolutions[i])
                Idents(object) = paste0("SCT_snn_res.",resolutions[i])
                UMAPPlot.1(object, group.by=paste0("SCT_snn_res.",resolutions[i]),pt.size = 0.3,label = T,
                           label.repel = T,alpha = 0.9,
                           do.return = F,
                           no.legend = T,label.size = 4, repel = T, 
                           title = paste("res =",resolutions[i],"in",g_name," based on harmony"),
                           do.print = T, save.path = save.path)
                file.rename(paste0(save.path,"UMAPPlot_object_SCT_snn_res.",resolutions[i],".jpeg"),
                            paste0(save.path,i,"-UMAPPlot_object_SCT_snn_res.",resolutions[i],".jpeg"))
                Progress(i,length(resolutions))
        }
}

# serial resolution and ICA heatmap on integrated data
if(step == 4){
        save.path <- paste0(path,g,"/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
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
                           title = paste("res =",res[i],"in",g_name," based on ICA"),
                           unique.name = "groups",
                           do.print = T, save.path = save.path)
                Progress(i,length(res))
        }
        
        object <- RunICA(object, verbose =F,nics = 100)
        a <- seq(1,50, by = 6)
        b <- a+5
        for(i in seq_along(a)){
                jpeg(paste0(save.path,"DimHeatmap_ica_",g,"_",a[i],"_",min(b[i],50),".jpeg"), units="in", width=10, height=7,res=600)
                DimHeatmap(object, dims = a[i]:min(b[i],50),
                           nfeatures = 30,reduction = "ica")
                dev.off() 
        }
}