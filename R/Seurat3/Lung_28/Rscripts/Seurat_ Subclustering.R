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
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
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

groups <- c("AT","B","En","F","Mon","SAE","SMG","SMP","T")
(g <- groups[args])
# === load data =============
(load(file = "data/Lung_28_20200116.Rda"))
object$cell_types %<>% gsub("-.*","",.) %>% gsub("[0-9]+","",.)
object$groups <- plyr::mapvalues(x = object$cell_types,
                                 from = c("AT","B","Car","D","DC",
                                          "En","F","Mac","MC","Mon",
                                          "NEC","Neu","Nr","P","PC",
                                          "Per","SAE","SM","SMG","T"),
                                 to = c("AT","B","Car","D","Mon",
                                        "En","F","Mon","MC","Mon",
                                        "SAE","Neu","Nr","P","PC",
                                        "SMP","SAE","SMP","SMG","T"))
Idents(object) = "groups"
object %<>% subset(idents = g)

#======1.6 without intergration =========================
DefaultAssay(object) = "SCT"
VariableFeatures(object) = rownames(object@assays$SCT@scale.data)
object <- RunICA(object, verbose =F,nics = 100)

npcs =100
object %<>% FindNeighbors(reduction = "ica",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.6)
object %<>% RunTSNE(reduction = "ica", dims = 1:npcs, check_duplicates = FALSE)
object %<>% RunUMAP(reduction = "ica", dims = 1:npcs)

TSNEPlot.1(object, group.by = "RNA_snn_res.0.6",label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F, unique.name = "groups",
           title = paste("Clusters in", g))
UMAPPlot.1(object, group.by = "RNA_snn_res.0.6",label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F, unique.name = "groups",
           title = paste("Clusters in", g))
save(object, file = paste0("data/Lung_28_",g,"_20200121.Rda"))

