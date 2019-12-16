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

Cell.types = c("labels","major_labels","family_labels","group_labels")
(cell.type <- Cell.types[args])

# load data
(load(file = "data/Lung_16_distal_20191022.Rda"))
object@meta.data$major_labels = gsub(":.*","",object@meta.data$labels)
object@meta.data$family_labels = object@meta.data$major_labels

object$family_labels[grep("Alveolar macrophages|Alveolar type 1|Alveolar type 2",
                          object$family_labels)] = "Alveolar epithelial cells"
object$family_labels[grep("Basal cells|Intermediate cells",
                          object$family_labels)] = "Basal cells+Intermediate cells"
object$family_labels[grep("Secretory cells|Mucus-producing|Serous cells",
                          object$family_labels)] = "Secretory family"
object$family_labels[grep("T cells|NK cells", object$family_labels)] = "T + NK cells"
object$family_labels[grep("Smooth muscle|Pericytes", object$family_labels)] = "Smooth muscle family"

object@meta.data$group_labels = object@meta.data$family_labels
object$group_labels[grep("Basal cells\\+Intermediate cells|Secretory family|Ciliated cells|Ionocytes/NEC|Neuroendocrine",
                         object$group_labels)] = "All airway epithelial cells"
object$group_labels[grep("Alveolar epithelial cells|Macrophages|Dendritic cells",
                          object$group_labels)] = "Macrophages and dendritic cells"

# subset
Idents(object) = "labels"
object <- subset(object, idents = c("Unknown","T cells:7SK.2+"),invert = T)
Idents(object) = cell.type
object %<>% sortIdent()
if(cell.type == "family_labels"){
        object <- subset(object, idents = c("Alveolar epithelial cells",
                                           "Basal cells+Intermediate cells",
                                           "Secretory family",
                                           "Smooth muscle family",
                                           "T + NK cells"))
}

if(cell.type == "group_labels"){
        object <- subset(object, idents = c("All airway epithelial cells",
                                                "Macrophages and dendritic cells"))
}

(labels = unique(object@meta.data[,cell.type]))
object@meta.data[,paste0(cell.type,"_conditions")] = paste0(object@meta.data[,cell.type],"_",
                                            as.character(object@meta.data$conditions))
Idents(object) = paste0(cell.type,"_conditions")
# Differential analysis
DefaultAssay(object) = "RNA"
Lung_markers <- FindPairMarkers(object,
                                ident.1 = paste0(labels,"_COPD"), 
                                ident.2 = paste0(labels,"_distal"),
                                logfc.threshold = 0.05, only.pos = F,
                                min.cells.group = 1,
                                test.use = "MAST",
                                save.files = FALSE)
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
write.csv(Lung_markers,paste0(path,"Lung_16_distal_COPD_",cell.type,"_conditions",".csv"))
Top_n = 5
top = Lung_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
object %<>% ScaleData(features=unique(c(as.character(top$gene))))

DoHeatmap.1(object, marker_df = Lung_markers, Top_n = 5, do.print=T,
            angle = 0,group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
            assay = "RNA",
            label=T, cex.row=5, legend.size = 5,width=10, height=7,unique.name = T,
            title = paste("Top 5 markers of",cell.type,"in distal and COPD sampels"))
