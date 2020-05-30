########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   
                   "magrittr","MAST","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
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

step = 1
# ==================================================
# generate ElbowPlot
if(step == 1){
        DEG_path <- paste0(path,g,"/DEG/")
        if(!dir.exists(DEG_path))dir.create(DEG_path, recursive = T)
        load(file = paste0("data/Lung_28_",g,"_20200126.Rda"))
        DefaultAssay(object) <- 'integrated'
        if(g == "AT") {
                object %<>% FindClusters(resolution = 0.03)
                object %<>% subset(idents = 0)
                res = 0.1
        }
        if(g == "B") {
                res = c(0.1)
                object %<>% FindClusters(resolution = res[1])
        }
        if(g == "D") {
                res = c(0.6)
                object %<>% FindClusters(resolution = res[1])
        }
        if(g == "En") {
                res = c(0.3)
                object %<>% FindClusters(resolution = res[1])
        }
        if(g == "F") {
                res = c(0.5,1.5)
                object %<>% FindClusters(resolution = res[2])
                Idents(object) = paste0("integrated_snn_res.",res[2])
                UMAPPlot.1(object, group.by = paste0("integrated_snn_res.",res[2]),
                           do.print = T, do.return = F,label = T, label.repel = T,
                           save.path = DEG_path)
                markers <- FindAllMarkers.UMI(object,assay = "SCT", logfc.threshold = 0.1,only.pos = T)
                write.csv(markers, file = paste0(DEG_path, "DEGs_",g,"_res.",res[2],".csv"))
                object %<>% FindClusters(resolution = res[1])
        }
        if(g == "Mon") {
                res = c(0.2)
                object %<>% FindClusters(resolution = res[1])
                object %<>% subset(idents = 7, invert = T)
        }
        if(g == "SAE") {
                res = c(0.3,0.8,1.5)
                for(i in 2:3){
                        object %<>% FindClusters(resolution = res[i])
                        Idents(object) = paste0("integrated_snn_res.",res[i])
                        UMAPPlot.1(object, group.by = paste0("integrated_snn_res.",res[i]),
                                   do.print = T, do.return = F,label = T, label.repel = T,
                                   save.path = DEG_path)
                        markers <- FindAllMarkers.UMI(object,assay = "SCT", logfc.threshold = 0.1,only.pos = T)
                        write.csv(markers, file = paste0(DEG_path, "DEGs_",g,"_res.",res[i],".csv"))
                }
                object %<>% FindClusters(resolution = res[1])
        }
        if(g == "SMG") {
                res = c(0.9)
                object %<>% FindClusters(resolution = res[1])
        }
        if(g == "SMP") {
                res = c(0.1,0.3)
                for(i in 1:2){
                        object %<>% FindClusters(resolution = res[i])
                }
                Idents(object) = paste0("integrated_snn_res.",res[2])
                UMAPPlot.1(object, group.by = paste0("integrated_snn_res.",res[2]),
                           do.print = T, do.return = F,label = T, label.repel = T,
                           save.path = DEG_path)
                markers <- FindAllMarkers.UMI(object,assay = "SCT", logfc.threshold = 0.1,only.pos = T)
                write.csv(markers, file = paste0(DEG_path, "DEGs_",g,"_res.",res[2],".csv"))
        }
        if(g == "T") {
                res = c(0.1,0.2)
                for(i in 1:2){
                        object %<>% FindClusters(resolution = res[i])
                }
                Idents(object) = paste0("integrated_snn_res.",res[2])
                UMAPPlot.1(object, group.by = paste0("integrated_snn_res.",res[2]),
                           do.print = T, do.return = F,label = T, label.repel = T,
                           save.path = DEG_path)
                markers <- FindAllMarkers.UMI(object,assay = "SCT", logfc.threshold = 0.1,only.pos = T)
                write.csv(markers, file = paste0(DEG_path, "DEGs_",g,"_res.",res[2],".csv"))
        }
        Idents(object) = paste0("integrated_snn_res.",res[1])
        UMAPPlot.1(object, do.print = T, do.return = F,label = T, label.repel = T,
                   save.path = DEG_path)
        DefaultAssay(object) = "SCT"
        markers <- FindAllMarkers.UMI(object,assay = "SCT", logfc.threshold = 0.1,only.pos = T)
        write.csv(markers, file = paste0(DEG_path, "DEGs_",g,"_res.",res[1],".csv"))
        #==== cells in each cluster ===== 
        df <- table(object@meta.data[,paste0("integrated_snn_res.",res[1])], object$orig.ident) %>% 
                as.data.frame()
        colnames(df) = c("cluster","samples","Freq")
        df %<>% spread("cluster","Freq")
        df$Region = gsub("-R","",df$samples) %>% gsub(".*-","",.)
        df = df[order(df$samples),]
        df = df[order(match(df$Region, c("P","D","T"))),]
        df %<>% select("Region", everything())
        write.csv(df, file = paste0(DEG_path, "Cell.num_",g,"_res.",res[1],".csv"))
        #==== genes in each cluster ===== 
        df <- table(object@meta.data[,paste0("integrated_snn_res.",res[1])], object$orig.ident) %>% 
                as.data.frame()
        colnames(df) =c("cluster","samples","nGene")
        
        df$cluster %<>% as.character()
        df$samples %<>% as.character()
        
        for(i in seq_len(nrow(df))) {
                cells <- object@meta.data[,paste0("integrated_snn_res.",res[1])] %in% 
                        df[i,"cluster"] & object$orig.ident %in% df[i,"samples"]
                df[i,"nGene"] = as.integer(mean(object@meta.data[cells,"nFeature_RNA"]))
        }
        df %<>% spread("cluster","nGene")
        df[is.na(df)] = 0
        df = df[order(df$samples),]
        df$Region = gsub("-R","",df$samples) %>% gsub(".*-","",.)
        df = df[order(match(df$Region, c("P","D","T"))),]
        df %<>% select("Region", everything())
        write.csv(df, file = paste0(DEG_path, "Gene.num_",g,"_res.",res[1],".csv"))
}
