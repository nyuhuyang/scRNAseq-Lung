########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#conda activate r4.0.3 
# devtools::install_github('satijalab/seurat-data') #3.1.5.9900
# remotes::install_github("mojaveazure/seurat-disk")
invisible(lapply(c("Seurat","dplyr","cowplot","magrittr","tidyr","tibble","data.table"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save_path <- "Yang/Lung_30/hg38/Deconvolution/"
if(!dir.exists(save_path))dir.create(save_path, recursive = T)

# 5.1.1 load data
# Rename ident
object <- readRDS(file = "data/Lung_SCT_30_20210831.rds") 
meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")
table(rownames(object@meta.data) == rownames(meta.data))
table(object$barcode ==meta.data$barcode)
object@meta.data = meta.data
object %<>% subset(subset =  Doublets %in% "Singlet" &
                           Cell_subtype != "Un")


tpm = data.table::fread(file = "data/RNA-seq/GTEx-Lung-tpm.txt",header = T)
genes = tpm[["GeneSymbol"]]

table(rownames(object) %in% genes)
table(genes %in% rownames(object))
share_genes = rownames(object)[rownames(object) %in% genes]
tpm %<>% filter(GeneSymbol %in% share_genes)
fwrite(tpm, file=paste0(save_path,"GTEx-Lung-tpm_hg38.txt"),quote = FALSE,
       row.names = FALSE,sep = "\t")

# generate reference
groups = c("Cell_subtype","Cell_type","UMAP_land","Family","Superfamily")
for(g in groups){
        print(g)
        Idents(object) = g
        expression = AggregateExpression(object,assay = "SCT",group.by = g)
        expression = expression$SCT
        tpm.mat = t(t(expression)*1e6 /colSums(expression) )
        write.table(tpm.mat[share_genes,],file = paste0(save_path,"psudobulk_tpm_Lung30_hg38_",g,"_reference.txt"),sep = "\t",quote = FALSE)
}


# add "gene" to the first line

# generate single-cell reference
counts = object[["SCT"]]@counts
cell_types <- sort(unique(object$cell_types))
meta.data = object@meta.data
select_barcode  <- pbapply::pbsapply(cell_types, function(cell_type){
        sub_meta.data = meta.data %>% filter(cell_types %in% cell_type)
        sample(rownames(sub_meta.data), size = min(nrow(sub_meta.data),100))
        
})
select_barcode %<>% unlist()
sub_object = subset(object,cells = select_barcode)

counts = as.matrix(sub_object[["SCT"]]@counts)
sub_meta.data =sub_object@meta.data
colnames(counts) = sub_object$cell_types
write.table(counts,file = paste0(path,"single.cell_Lung30_cell.type_reference.txt"),sep = "\t",quote = FALSE)

# add "GeneSymbol" to the first line.

# combine CIBERSORT with meta.data
GTE_meta.data = read.csv("data/RNA-seq/GTEx Portal - lung sample info.csv",
                         stringsAsFactors = F)


CIBERSORTx_GTEx_list <- list.files(save_path, pattern = "GTEx-Lung-tpm_by_Lung30_hg38_.*xlsx",full.names = T)
CIBERSORTx_GTEx <- pbapply::pblapply(CIBERSORTx_GTEx_list, function(x){
        tmp <- readxl::read_excel(x)
        colnames(tmp) %<>% sub("Mixture", "Tissue.Sample.ID", .)
        tmp$Tissue.Sample.ID %<>% gsub("-SM.*","",.)
        table(tmp$Tissue.Sample.ID %in% GTE_meta.data$Tissue.Sample.ID)
        tmp = left_join(tmp, GTE_meta.data, by = "Tissue.Sample.ID")
        tmp
})
GTEx_list = CIBERSORTx_GTEx_list %>% gsub(".*GTEx-Lung-tpm_by_Lung30_hg38_","",.) %>% gsub("_Results.*","",.)
names(CIBERSORTx_GTEx) = GTEx_list
#CIBERSORTx_GTEx = CIBERSORTx_GTEx[c("super.family","cell.family","RS.group","major.cell_types")]
openxlsx::write.xlsx(CIBERSORTx_GTEx, file =  paste0(save_path,"CIBERSORTx-GTEx_summary.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#
