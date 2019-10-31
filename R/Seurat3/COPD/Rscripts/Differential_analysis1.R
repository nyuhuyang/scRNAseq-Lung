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

conditions = c("COPD","distal","conditions","RNA_snn_res.0.8")
(con <- conditions[args])

# load data
(load(file = "data/Lung_16_distal_20191022.Rda"))

if(con %in% c("COPD","distal")){
        Idents(object) = "conditions"
        object <- subset(object, idents = con)
        Idents(object) = "labels"
} else Idents(object) = con

object %<>% sortIdent(numeric = T)
table(Idents(object))
Lung_markers <- FindAllMarkers.UMI(object,logfc.threshold = 0.5,only.pos = T, 
                                          min.pct = 0.1,return.thresh = 0.05)
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
write.csv(Lung_markers,paste0(path,"Lung_16_",con,"_markers.csv"))
Top_n = 5
top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
object %<>% ScaleData(features=unique(c(as.character(top$gene))))
featuresNum <- make.unique(features, sep = ".")
object %<>% MakeUniqueGenes(features = features)

DoHeatmap.1(object, features = featuresNum, Top_n = Top_n, do.print=T,
            angle = 0,group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
            assay = "RNA",
            label=T, cex.row= 0.1, legend.size = 5,width=10, height=7,unique.name = T,
            title = paste("Top",Top_n,"markers of",con,"in 16 D+COPD sampels"))
