########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   
                   "magrittr","MAST","future","ggplot2","tidyr"), function(x) {
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

sheets = c("T20","IE", "DS", "3C")
(sheet <- groups[args])
df_samples <- readxl::read_excel("doc/Cell type markers for UMAP re-clustering.xlsx",
                                 sheet = sheet)

save.path <- paste0(path,sheet,"/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)
# ==================================================
step = 1
# Find pc number
if(step == 1){
        #======1.2 load  Seurat =========================
        (load(file = "data/Lung_28_20200102.Rda"))
        object[["integrated"]] = NULL
        object@neighbors = list()
        object@reductions = list()
        DefaultAssay(object) = "SCT"
        object <- FindVariableFeatures(object = object, selection.method = "vst",
                                       num.bin = 20,nfeatures = 10000,
                                       mean.cutoff = c(0.1, 8), 
                                       dispersion.cutoff = c(1, Inf))
        # keep the VariableFeatures in marker list only
        marker_genes = FilterGenes(object, df_samples$genes)
        VariableFeatures(object) %<>% .[. %in% marker_genes]
        length(VariableFeatures(object))
        # Identify the 20 most highly variable genes
        top20 <- head(VariableFeatures(object), 20)
        # plot variable features with and without labels
        plot1 <- VariableFeaturePlot(object)
        plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
        jpeg(paste0(save.path,"VariableFeaturePlot.jpeg"), units="in", width=10, height=7,res=600)
        print(plot2)
        dev.off()
        hvf.info <- HVFInfo(object = object)
        hvf.info = hvf.info[VariableFeatures(object),]
        write.csv(hvf.info, file = paste0(save.path,"high_variable_genes.csv"))
        
        #======1.7 harmony =========================
        object %<>% ScaleData
        object %<>% RunPCA(verbose = T,npcs = 100, features = VariableFeatures(object))
        PCAPlot.1(object, do.print = T, save.path = save.path)
        
        npcs= 100
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
        saveRDS(object, file = paste0("data/Lung_28_",args,"-",sheet,"-harmony_2020330.rds"))
}