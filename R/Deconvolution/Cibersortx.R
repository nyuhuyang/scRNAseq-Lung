########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#conda activate r4.0.3 
# devtools::install_github('satijalab/seurat-data') #3.1.5.9900
# remotes::install_github("mojaveazure/seurat-disk")
invisible(lapply(c("Seurat","dplyr","cowplot","magrittr","tidyr","tibble"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save_path <- "Yang/Lung_30/hg19/Deconvolution/"
if(!dir.exists(save_path))dir.create(save_path, recursive = T)

# 5.1.1 load data
# Rename ident
object = readRDS(file = "data/Lung_SCT_30_20200710.rds")
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
Idents(object) = "cell_types"
object %<>% subset(idents = "Un", invert = T)
df_samples <- readxl::read_excel("doc/Annotations/Cell type abbreviation - refs for deconvolution.xlsx")
table(df_samples$cell_types %in% object@meta.data[,"cell_types"])

groups = c("cell_types.select","cell_types","major.cell_types","RS.group","cell.family","super.family")
for(g in groups){
        object[[g]] = plyr::mapvalues(object@meta.data[,"cell_types"],
                                      from = df_samples$cell_types,
                                      to = df_samples[[g]])
}

for(g in groups){
        freq_table <- prop.table(x = table(object@meta.data[[g]],object$orig.ident),
                                 margin = 2) %>% as.data.frame.matrix
        write.table(t(freq_table), paste0(save_path,"psudobulk_",g,"_percent.txt"),sep = "\t",quote = FALSE)
}
# remove "---" column
# generate mixture
Idents(object) = "orig.ident"
expression = AggregateExpression(object,assay = "SCT")
expression = expression$SCT
tpm.mat = t(t(expression)*1e4 /colSums(expression) )
write.table(tpm.mat,file = paste0(save_path,"psudobulk_tpm_Lung30_orig.ident.txt"),sep = "\t",quote = FALSE)
# add "GeneSymbol" to the first line.

# convert to h5ad
object[["pca"]] = NULL
SaveH5Seurat(object, filename = "data/lung_30.h5Seurat")
Convert("data/lung_30.h5Seurat", dest = "h5ad")


meta_color = table(object$cell_types, object$cell_types.colors) %>% t %>% as.data.frame.matrix()
((meta_color > 0)*1) %>% rowSums()
UMAPPlot(object,reduction = NULL, cols  =ExtractMetaColor(object)) + NoLegend()

meta.data = object@meta.data
meta.data = meta.data[!duplicated(meta.data$cell_types.colors),]

freq_table = freq_table *100
Cell_proportion <- readxl::read_excel(paste0("doc/","Cell proportion - single cell data.xlsx"))
table(Cell_proportion$cell.types %in% rownames(freq_table))
table(colnames(Cell_proportion) %in% colnames(freq_table))

freq_table = freq_table[Cell_proportion$cell.types,colnames(Cell_proportion)[-1]]
write.table(freq_table, paste0(path,"psudobulk_cell_type_percentage.txt"),sep = "\t",quote = FALSE)
freq_table$cell_types = rownames(freq_table)

freq_table %<>% gather("sample","percentage",-cell_types)
freq_table$cell_types.color = plyr::mapvalues(freq_table$cell_types, from  = meta.data$cell_types,
                                              to = meta.data$cell_types.colors)
ggOut = ggplot(freq_table, aes_string(x = "sample",
                                      y =    "percentage",
                                      fill = "cell_types")) +
        geom_col() + 
        ggtitle("ground truth")+
        xlab("sample")+
        ylab("Cell Proportion (%)")+#NoLegend()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size =10),
              plot.title = element_text(hjust = 0.5))+
        scale_fill_manual(values = ExtractMetaColor(object,group.by = "cell_types"))
ggOut

cross.entropy <- function(p, phat){
        x <- 0
        for (i in 1:length(p)){
                if(phat[i] == 0) next
                x <- x + (p[i] * log(phat[i]))
        }
        return(-x)
}
cross.entropy(freq_table[freq_table$sample %in% "CU-12-D","percentage"],
              freq_table[freq_table$sample %in% "CU-12-D","percentage"]
)


# generate reference
for(g in groups){
        print(g)
        Idents(object) = g
        if("---" %in% unique(Idents(object))){
                sub_object = subset(object, idents = "---", invert = T)
        } else sub_object = object
        
        expression = AggregateExpression(sub_object,assay = "SCT",group.by = g)
        expression = expression$SCT
        tpm.mat = t(t(expression)*1e6 /colSums(expression) )
        write.table(tpm.mat,file = paste0(save_path,"psudobulk_tpm_Lung30_",g,"_reference.txt"),sep = "\t",quote = FALSE)
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


CIBERSORTx_GTEx_list <- list.files(save_path, pattern = "CIBERSORTx_GTEx.*_Results.xlsx",full.names = T)
CIBERSORTx_GTEx <- pbapply::pblapply(CIBERSORTx_GTEx_list, function(x){
        tmp <- readxl::read_excel(x)
        colnames(tmp) %<>% sub("Mixture", "Tissue.Sample.ID", .)
        tmp$Tissue.Sample.ID %<>% gsub("-SM.*","",.)
        table(tmp$Tissue.Sample.ID %in% GTE_meta.data$Tissue.Sample.ID)
        tmp = left_join(tmp, GTE_meta.data, by = "Tissue.Sample.ID")
        tmp = tmp[order(tmp$Sex,tmp$Age.Bracket),]
})
GTEx_list = CIBERSORTx_GTEx_list %>% gsub(".*CIBERSORTx_GTEx-","",.) %>% gsub("_Results.*","",.)
names(CIBERSORTx_GTEx) = GTEx_list
CIBERSORTx_GTEx = CIBERSORTx_GTEx[c("super.family","cell.family","RS.group","major.cell_types","cell_types","selected.cell_types")]
openxlsx::write.xlsx(CIBERSORTx_GTEx, file =  paste0(save_path,"/CIBERSORTx-GTEx_summary.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#
