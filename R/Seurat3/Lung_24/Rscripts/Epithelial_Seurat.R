########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","magrittr",
                   "tidyr","gplots","MAST"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
set.seed=101
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# conditions
conditions = c("proximal","distal","terminal","All")
(con <- conditions[args])

# Load Seurat
(load(file="data/Lung_24_20191206.Rda"))
table(object$orig.ident)
object$major_cell.types = gsub(":.*","",object$cell.types)
Idents(object) = "major_cell.types"
object %<>% subset(idents = c("Alveolar type 1 cells","Alveolar type 2 cells","Basal cells",
                              "Ciliated cells","Hybrid cells","Intermediate cells",
                              "Ionocytes","Mucus-producing cells","Myoepithelial cells",
                              "Neuroendocrine cells","Pre-ciliated cells","Secretory cells",
                              "Squamous","Submucosal gland"))

Idents(object) = "cell.types"
if (con == "proximal") {
    object %<>% subset(idents = c("Alveolar type 1 cells","Alveolar type 2 cells:A",
                                  "Alveolar type 2 cells:B","Alveolar type 2 cells:C",
                                  "Secretory cells:Distal","Secretory cells:Distal:2"), 
                       invert = TRUE)
}
if (con == "distal") {
    object %<>% subset(idents = c("Myoepithelial cells","Squamous",
                                  "Submucosal gland:Mucous cells",
                                  "Submucosal gland:Serous cells"), 
                       invert = TRUE)
}
if (con == "terminal") {
    object %<>% subset(idents = c("Mucus-producing cells","Myoepithelial cells",
                                  "Squamous","Submucosal gland:Mucous cells",
                                  "Submucosal gland:Serous cells"),
                       invert = TRUE)
}
Idents(object) = "conditions"
if(con != "All") object <- subset(object, idents = con)
######################################
DefaultAssay(object) = "SCT"
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(object), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
jpeg(paste0(path,"VariableFeaturePlot_",con,".jpeg"), units="in", width=10, height=7,res=600)
print(plot2)
dev.off()
#======1.3 1st run of pca-tsne  =========================
npcs =50
object <- RunPCA(object, features = VariableFeatures(object),verbose =F,npcs = npcs)

object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
#object %<>% RunTSNE(reduction = "pca", dims = 1:npcs, check_duplicates = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

object@assays$RNA@scale.data = matrix(0,0,0)
save(object, file = paste0("data/Epi_24_",con,"_20191223.Rda"))
Idents(object) = "RNA_snn_res.0.8"
object %<>% sortIdent(numeric = T)
UMAPPlot.1(object, group.by="RNA_snn_res.0.8",pt.size = 1,label = T,
           label.repel = T,alpha = 0.9,
           no.legend = F,label.size = 4, repel = T, title = paste(con,"clusters with No Integration"),
           do.return = F,do.print = T,unique.name = "conditions")
Idents(object) = "cell.types"
object %<>% sortIdent
UMAPPlot.1(object, group.by="cell.types",pt.size = 1,label = T,
           label.repel = T,alpha = 0.9,cols = ExtractMetaColor(object),
           no.legend = F,label.size = 4, repel = T, title = paste(con,"cell types"),
           do.return = F,do.print = T,unique.name = "conditions")

