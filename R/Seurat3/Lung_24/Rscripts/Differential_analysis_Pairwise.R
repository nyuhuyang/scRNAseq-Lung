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

Pairs = c("P.vs.D","P.vs.(D+T combined)","D.vs.T")
labels = c("cell.types","major.cell.types","family.cell.types")
Pairs=paste0(labels,"_", rep(Pairs,each =3))
(p <- Pairs[args])

# load data
(load(file = "data/Lung_24_20190918.Rda"))
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
object@meta.data$cell.types %<>% gsub("Dendritic Cells","Dendritic cells",.)
# group cell types
object@meta.data$major.cell.types = gsub(":.*","",object@meta.data$cell.types)
object@meta.data$major.cell.types = gsub("-[0-9]","",object@meta.data$major.cell.types)

object$major.cell.types[grep("Mucus-producing cells",object$major.cell.types)] = "Secretory cells"
object$major.cell.types[grep("Pre-ciliated cells",object$major.cell.types)] = "Ciliated cells" 
object$family.cell.types = object@meta.data$major.cell.types

object$family.cell.types[grep("Basal.*|Secretory.*|Ciliated.*|Intermediate.*|Ionocytes|Neuroendocrine.*",
                          object$family.cell.types)] = "Surface Airway Epithelium"
object$major.cell.types[grep("Smooth.*|Pericytes",object$major.cell.types)] = "Smooth muscle cells & Pericytes" 
object$family.cell.types[grep("Fibroblasts|Smooth.*",
                          object$family.cell.types)] = "Stromal cells"

# subset
label = sub("\\_.*","",p)
pair = sub(".*\\_","",p)

Idents(object) = "conditions"
if(pair == "P.vs.D"){
        sub_object <- subset(object, idents = c("distal","proximal"))
}
if(pair == "P.vs.(D+T combined)"){
        sub_object <- object
        sub_object@meta.data$conditions %<>% gsub("distal|terminal","distal+terminal",.)
}
if(pair == "D.vs.T"){
        sub_object <- subset(object, idents = c("distal","terminal"))
}
if(label == "family.cell.types"){
        Idents(sub_object) = "family.cell.types"
        sub_object <- subset(sub_object, idents = c("Surface Airway Epithelium",
                                                "Stromal cells"))
}

cell.type = unique(sub_object@meta.data[,label]) %>% sort
sub_object@meta.data[,paste0(label,"_conditions")] = paste0(sub_object@meta.data[,label],"_",
                                            as.character(sub_object@meta.data$conditions))
Idents(sub_object) = paste0(label,"_conditions")

# remove cluster with less than 3 cells======
label_conditions <- as.factor(sub_object@meta.data[,paste0(label,"_conditions")])
levels = paste0(cell.type,"_",rep(unique(sub_object$conditions),
                                  each =length(cell.type)))
label_conditions %<>% factor(levels = levels)
table_subset <- table(label_conditions) %>% as.data.frame
(remove <- table_subset[table_subset$Freq < 3,"label_conditions"] %>% 
                gsub("\\_.*","",.) %>% unique() %>% as.character())
(cell.type <- cell.type[!(cell.type %in% remove)])
Idents(sub_object) = label
sub_object %<>% subset(idents= remove, invert = T)

# Differential analysis
DefaultAssay(object) = "RNA"
Idents(sub_object) = paste0(label,"_conditions")
sub_object %<>% sortIdent
if(pair == "P.vs.D"){
        Lung_markers <- FindPairMarkers(sub_object,
                                        ident.1 = paste0(cell.type,"_distal"), 
                                        ident.2 = paste0(cell.type,"_proximal"),
                                        logfc.threshold = 0.05, only.pos = F,
                                        test.use = "MAST",
                                        save.files = FALSE)
}
if(pair == "P.vs.(D+T combined)"){
        Lung_markers <- FindPairMarkers(sub_object,
                                        ident.1 = paste0(cell.type,"_distal+terminal"), 
                                        ident.2 = paste0(cell.type,"_proximal"),
                                        logfc.threshold = 0.05, only.pos = F,
                                        test.use = "MAST",
                                        save.files = FALSE)
}
if(pair == "D.vs.T"){
        Lung_markers <- FindPairMarkers(sub_object,
                                        ident.1 = paste0(cell.type,"_terminal"), 
                                        ident.2 = paste0(cell.type,"_distal"),
                                        logfc.threshold = 0.05, only.pos = T,
                                        test.use = "MAST",
                                        save.files = FALSE)
}
write.csv(Lung_markers,paste0(path,"Lung_24_",p,"~.csv"))
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
write.csv(Lung_markers,paste0(path,"Lung_24_",p,".csv"))
