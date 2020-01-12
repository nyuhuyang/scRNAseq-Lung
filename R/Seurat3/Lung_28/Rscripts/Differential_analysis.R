########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr",
                   "MAST","future","gplots"), function(x) {
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

if(args == 1){
        (load(file = "data/Lung_28_20200103.Rda"))
        DefaultAssay(object)  = "SCT"
        (load(file = paste0("output/SCINA_20200110.Rda")))
        object@meta.data$SCINA = results$cell_labels
        df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
        object$SCINA %<>% plyr::mapvalues(
                from = df_cell_types$`Cell types`,
                to = df_cell_types$Abbreviation)
        Idents(object) = "Doublets"
        object %<>% subset(idents = "Singlet")
        Idents(object) = "SCINA"
        object %<>% sortIdent()
        system.time(Lung_markers <- FindAllMarkers.UMI(object, logfc.threshold = 0.25, only.pos = T,
                                           test.use = "MAST"))
        Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
        write.csv(Lung_markers,paste0(path,"Lung_28-FC0.25_markers_SCINA.csv"))
}
