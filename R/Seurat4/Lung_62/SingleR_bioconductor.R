#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)
library(dplyr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# ====== load old Lung data=============
Lung = readRDS(file = "data/Lung_SCT_30_20200710.rds")

# ====== load new Lung data=============
object = readRDS(file = "data/Lung_SCT_59_20210814.rds")

Lung$barcode = gsub(".*_","",colnames(Lung)) %>% gsub("-1$","",.)

Lung$orig.ident = gsub("-","_",Lung$orig.ident)
Lung$orig.ident %<>% gsub("VU19_D","VU_19_D",.)
table(Lung$orig.ident %in% object$orig.ident)

Lung$barcode = paste0(Lung$orig.ident,"-",Lung$barcode)
table(Lung$barcode %in% colnames(object))

meta_data = Lung@meta.data[,c("barcode","cell_types","Doublets")]
rownames(meta_data) = meta_data$barcode
meta_data = meta_data[Lung$barcode %in% colnames(object),]
object$barcode = colnames(object)
table(meta_data$barcode %in% object$barcode)
meta.data = object@meta.data
meta.data$cell_types = NULL
meta.data %<>% left_join(meta_data,by = "barcode")
meta.data[is.na(meta.data$cell_types),"cell_types"] = "Unknown"

rownames(meta.data) = meta.data$barcode
table(rownames(meta.data) == rownames(object@meta.data))
object@meta.data = meta.data
sub_object = subset(object, subset = cell_types != "Un")
sub_object %<>% subset(subset = Doublets == "Singlet")
rm(Lung);GC()
#======== generate reference =================
df_samples <- readxl::read_excel("doc/Annotations/Cell type abbreviation - refs for deconvolution.xlsx")
table(df_samples$cell_types %in% sub_object@meta.data[,"cell_types"])
save_path = "Yang/Lung_59/Deconvolution/"
groups = c("cell_types.select","cell_types","major.cell_types","RS.group","cell.family","super.family")
for(g in groups){
    sub_object[[g]] = plyr::mapvalues(sub_object@meta.data[,"cell_types"],
                                  from = df_samples$cell_types,
                                  to = df_samples[[g]])
}

# generate reference
for(g in groups){
    print(g)
    Idents(sub_object) = g
    if("---" %in% unique(Idents(sub_object))){
        sub_object1 = subset(sub_object, idents = "---", invert = T)
    } else sub_object1 = sub_object
    
    expression = AggregateExpression(sub_object1,assay = "SCT",group.by = g)
    expression = expression$SCT
    tpm.mat = t(t(expression)*1e6 /colSums(expression) )
    write.table(tpm.mat,file = paste0(save_path,"psudobulk_tpm_Lung30_hg38_",g,"_reference.txt"),sep = "\t",quote = FALSE)
}

rm(sub_object);rm(sub_object1);GC()



sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# 2. check and prepare reference ==============================
Lung_bulk <- read.delim(file="Yang/Lung_59/Deconvolution/psudobulk_tpm_Lung30_cell_types.select_reference.txt")

meta.data = data.frame("label.fine" = colnames(Lung_bulk),
                       "label.ont" = colnames(Lung_bulk))
Lung_sce <- SummarizedExperiment(list(logcounts=Lung_bulk),
                                  colData=DataFrame(meta.data))

common <- Reduce(intersect, list(rownames(sce),
                                 rownames(Lung_sce)
))

table(Lung_sce$label.fine)
system.time(trained <- trainSingleR(ref = Lung_sce,
                                    labels=Lung_sce$label.fine))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 542.727 sec
saveRDS(object = pred, file = "output/Lung_SCT_59_20210814_singleR_pred.rds")
