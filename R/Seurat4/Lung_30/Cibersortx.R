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


tpm = data.table::fread(file = "data/RNA-seq/GTEx-Lung-tpm~.csv",header = T)
genes = tpm[["V1"]]

table(rownames(object) %in% genes)
table(genes %in% rownames(object))
share_genes = rownames(object)[rownames(object) %in% genes]
tpm %<>% filter(GeneSymbol %in% share_genes)
fwrite(tpm, file=paste0(save_path,"GTEx-Lung-tpm_hg38.txt"),quote = FALSE,
       row.names = FALSE,sep = "\t")

# generate raw single-cell reference
# read DE genes
csv_list <- list.files("output/20211022/ASE/",pattern = "_ASE.csv",full.names = T)
deg <- pbapply::pblapply(csv_list,function(csv){
        tmp = read.csv(csv,row.names = 1)
        tmp$gene = rownames(tmp)
        tmp$cluster = sub("_vs_ASE.csv","",csv) %>% sub(".*[0-9]+-","",.)
        tmp
}) %>% bind_rows()
#deg %<>%  group_by(cluster) %>% top_n(400, avg_log2FC)
length(unique(deg$gene))
table(deg$cluster)
fwrite(deg, file="output/20211022/ASE/ASE_DEGs.txt",quote = FALSE,
       row.names = FALSE,sep = "\t")
table(deg[deg$gene %in% deg[deg$cluster == "TASC","gene"],"cluster"])
table(deg[!deg$gene %in% deg[deg$cluster != "TASC","gene"],"cluster"])

Idents(object) = "Cell_subtype"
ASE = subset(object, idents = c("BC", "IC", "S1", "S-Muc", "TASC", "H", "p-C", "C1", "C-s", "Ion", "NE"))
mt = as.matrix(ASE[["SCT"]]@counts[unique(deg$gene),])
#keep_genes = rowSums(ASE[["SCT"]]@counts) ;table(keep_genes)
#keep_genes = keep_genes > 100
#mt = as.matrix(ASE[["SCT"]]@counts[keep_genes,])
dt = as.data.table(mt)
colnames(dt) = make.unique(ASE$Cell_subtype)
GeneSymbol = data.frame("GeneSymbol" =  rownames(mt))
dt = cbind(GeneSymbol,dt)

fwrite(dt, file=paste0(save_path,"ASE_scRNA_reference_hg38_full.txt"),quote = FALSE,
       row.names = FALSE,sep = "\t")

# generate psedo bulk reference
groups = c("Cell_subtype","Cell_type","UMAP_land","Family","Superfamily")
for(g in groups){
        print(g)
        Idents(object) = g
        expression = AggregateExpression(object,assay = "SCT",group.by = g)
        expression = expression$SCT
        tpm.mat = t(t(expression)*1e6 /colSums(expression) )
        GeneSymbol = data.frame("GeneSymbol" =  rownames(tpm.mat))
        tpm.mat = cbind(GeneSymbol,tpm.mat)
        rownames(tpm.mat) = rownames(tpm.mat)
        
        deg <- readxl::read_excel("Yang/Lung_30/hg38/DE_analysis/Lung_30_DEG_Cell.category.xlsx",sheet = g)
        keep_genes = unique(deg$gene)
        print(length(keep_genes))
        
        fwrite(tpm.mat[keep_genes,],file = paste0(save_path,"psudobulk_tpm_Lung30_hg38_",g,"_reference.txt"),quote = FALSE,
               row.names = FALSE,sep = "\t")
}
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
