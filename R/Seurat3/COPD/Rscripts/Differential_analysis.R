########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","tidyr","magrittr","gplots","MAST"), function(x) {
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

Cell.types = c("labels","major_labels")
(cell.type <- Cell.types[args])

# load data
(load(file = "data/Lung_16_distal_20191022.Rda"))

# Differential analysis
(labels = unique(object@meta.data[,cell.type]))
object@meta.data[,paste0(cell.type,"_conditions")] = paste0(object@meta.data[,cell.type],"_",
                                            as.character(object@meta.data$conditions))
Idents(object) = paste0(cell.type,"_conditions")
DefaultAssay(object) = "RNA"
Lung_markers <- FindPairMarkers(object,
                                ident.1 = paste0(labels,"_COPD"), 
                                ident.2 = paste0(labels,"_distal"),
                                logfc.threshold = 0.05, only.pos = F,
                                test.use = "MAST")
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
write.csv(Lung_markers,paste0(path,"Lung_16_distal_COPD_markers.csv"))
Top_n = 5
top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
object %<>% ScaleData(features=unique(c(as.character(top$gene))))

DoHeatmap.1(object, marker_df = Lung_markers, Top_n = 5, do.print=T,
            angle = 0,group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
            assay = "RNA",
            label=T, cex.row=5, legend.size = 5,width=10, height=7,unique.name = T,
            title = paste("Top 5 markers of",cell.type,"in distal and COPD sampels"))