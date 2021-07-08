#conda activate r4.0
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","magrittr","caret","HDF5Array",
                   "Matrix","hdf5r"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
save.path = "Yang/Lung_30/Feature_selection/"
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

########################################################################
#
#  1. load data and prepare training set
# 
# ######################################################################

object = readRDS(file = "data/Lung_30_20200710.rds") 
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
object$cell_types %<>% gsub("d-S","TASC",.)
SAE <- c("BC1","BC2","BC-p","IC1","IC2","IC3","S","TASC","H","p-C","C1","C2","C3","Ion","NE")
Idents(object) = "cell_types"
object %<>% subset(idents = SAE)
Idents(object) = "conditions"

#Step 1: Establish region-specific signatures
# Compared groups: P vs D (enhanced option P vs D+T)
object %<>% subset(idents =c("proximal","distal"))
#Split object meta.data into 2 group =====
split_meta.data <- split(object@meta.data, f = object$conditions)
#Split each meta.data into k group =====
lapply(split_meta.data,function(x) length(x$conditions))
k = 1000 
ident = "conditions"
continuous.label = c("proximal","distal")
for(i in 1:length(split_meta.data)){
        meta_index <- caret::createFolds(split_meta.data[[i]][,ident],
                                         k = k, list = TRUE, 
                                         returnTrain = FALSE)
        for(n in 1:length(meta_index)){
                split_meta.data[[i]][meta_index[[n]],"caret"] = 
                        paste(split_meta.data[[i]][meta_index[[n]],ident],
                              n, sep = "_")
        }
}
meta.data = bind_rows(split_meta.data)
object@meta.data = meta.data[rownames(object@meta.data),]
object$caret %<>% factor(levels = paste0(rep(continuous.label,each = k),
                                           "_",rep(1:k)))
Idents(object) = "caret"
expr <- AverageExpression(object, assays = "SCT")
dataset = expr[[1]]
table(rowSums(dataset) > 10)

dataset = dataset[rowSums(dataset) > 10,]
y = gsub("_.*","",colnames(dataset)) %>% plyr::mapvalues(from = c("proximal","distal"),
                                                         to = c(0,1)) %>% as.numeric()
dataset = rbind(dataset,Class = y)
dim(dataset)
dataset[(nrow(dataset)-4):nrow(dataset),c(1:2,1999:2000)]
#dataset <- Matrix(dataset, sparse = TRUE)
# save mat
h5file = paste0(save.path,"P_D_mat.h5")
sm <- writeHDF5Array(dataset, h5file, name="data")
sm[(nrow(dataset)-4):nrow(dataset),c(1:2,1999:2000)]
chunkdim(sm)
h5ls(h5file)
h5writeDimnames(dimnames(dataset), h5file, "data")
h5ls(h5file)
h5closeAll()


#=============
featureScores = read.csv("Yang/Lung_30/Feature_selection/chi2.csv",)
featureScores = featureScores[order(featureScores$Score,decreasing = T),]
head(featureScores)
