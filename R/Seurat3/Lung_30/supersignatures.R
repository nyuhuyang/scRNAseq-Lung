invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots","data.table"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
#========
read.path = "Yang/Lung_30/DE_analysis/"
list_files <- list.files(path = paste0(read.path,"C_Cell_types"),
                         pattern ="Lung_30-")
int <- gsub("Lung_30-","",list_files) %>% gsub("_.*","",.) %>% as.integer()
table(1:62 %in% int)

list_files_C <- list.files(path = paste0(read.path,"C_Cell_types"),
                           pattern ="Lung_30-",full.names = T)

deg_list = pbapply::pblapply(list_files_C, function(x) {
         read.csv(x,row.names = 1, stringsAsFactors = F) %>% 
                mutate(pct.1_pct.2 = pct.1-pct.2) %>%
                filter(avg_logFC >=1) %>%
                filter(pct.1 >=0.5) %>%
                #filter(pct.2 < 0.3) %>%
                filter(pct.1_pct.2 > 0.4)

})
supersignatures = bind_rows(deg_list)
rownames(supersignatures) = make.unique(supersignatures$gene)
anno <- readxl::read_excel("doc/Annotations/Cell type abbreviation.xlsx")
supersignatures$cluster %<>% plyr::mapvalues(from = anno$Abbreviation,
                                             to = anno$`Revised abbreviations`)
write.csv(supersignatures, 
          file =  paste0(paste0(read.path,"C_Cell_types/supersignatures_Extended.csv")))


# combine EVG
read.path = "Yang/Lung_30/DE_analysis/F_EVGs_allCells/"
superfamily <- c("Epithelial","Structural","Immune")

EVGs_list <- lapply(superfamily,function(s) {
        readxl::read_excel(paste0(read.path,"Lung_30-EVGs.xlsx"), sheet = s)
})
cbind.fill <- function(...){
        nm <- list(...) 
        nm <- lapply(nm, as.matrix)
        n <- max(sapply(nm, nrow)) 
        do.call(cbind, lapply(nm, function (x) 
                rbind(x, matrix("", n-nrow(x), ncol(x))))) 
}
EVGs_df <- do.call(cbind.fill, EVGs_list)
EVGs_df = apply(EVGs_df,2,as.character)

deg_list <- split(supersignatures, f = supersignatures$cluster)
cell.types = names(deg_list)
deg_list1 = pbapply::pblapply(deg_list, function(x) {
                          cell.type = unique((x$cluster))
                          x$EVG = x$gene %in% na.omit(EVGs_df[,cell.type])
                          x
                 })

deg_list1 <- list()
for(i in seq_along(deg_list)) {
        x = deg_list[[i]]
        cell.type = unique((x[,"cluster"]))
        if(!(cell.type %in% colnames(EVGs_df))) {
                x[,"EVG"] = FALSE
                next
                }
        x[,"EVG"] = x[,"gene"] %in% na.omit(EVGs_df[,cell.type])
        deg_list1[[i]] = x
}
supersignatures1 = bind_rows(deg_list1)
supersignatures1 = supersignatures1[order(supersignatures1$gene),]
supersignatures1$EVG  %<>% plyr::mapvalues(from = c(TRUE, FALSE),
                                           to = c("yes","no"))
read.path = "Yang/Lung_30/DE_analysis/"
write.csv(supersignatures1, 
          file =  paste0(paste0(read.path,"C_Cell_types/intersection_supersignatures_Extended.csv")))

# =======

df_TPM = fread("data/RNA-seq/GTEx-Lung-tpm~.csv",header = T)
TPM <- column_to_rownames(df_TPM, var = "V1")
dim(TPM)
# or 
load("data/Lung_GTEx_20210226.Rda")
TPM <- as.matrix(object[["RNA"]]@counts)
dim(TPM)

temp = readxl::read_excel("Yang/GTEx/Cell type EVG genes - 577 GTEx samples ordered.xlsx")
temp = temp[-1]

# replace gene name 
#all_genes = unique(supersignatures[,"gene"])
#old_gene = all_genes[!(all_genes %in% rownames(TPM))]
#"PLAAT4" %in% rownames(TPM)


exp_list <- list()
for(i in seq_along(cell.types)) {
        genes <- supersignatures[supersignatures$cluster == cell.types[i],"gene"]
        genes = genes[genes != ""]
        genes = genes[!is.na(genes)]
        genes = genes[(genes %in% rownames(TPM))]
        exp <- TPM[genes,colnames(temp)]
        exp_list[[i]] = exp#rbind.data.frame(temp,exp)
        colnames(exp_list[[i]])[1] = ""
        Progress(i, length(cell.types))
}

names(exp_list) = cell.types
#openxlsx::write.xlsx(exp_list, file =  "Yang/GTEx/Cell type supersignatures genes - 577 GTEx samples ordered.xlsx",
#                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))


# final integrated signatures=================
Int_supersignatures <- readxl::read_excel("doc/Integrated-signatures-categorized.xlsx") 
sheets = unique(Int_supersignatures$`Excel sheet`)

# replace gene name 
all_genes = pull(Int_supersignatures[,"Gene"])
table(all_genes %in% rownames(TPM))
old_gene = all_genes[!(all_genes %in% rownames(TPM))]
write.csv(old_gene,file = paste0(path, "old_genes.csv"))

Int_supersignatures %<>% as.data.frame()
exp_list <- list()
for(i in seq_along(sheets)) {
        row_tmp = Int_supersignatures[Int_supersignatures$`Excel sheet` == sheets[i],]
        genes <- row_tmp$Gene
        genes = genes[genes != ""]
        genes = genes[!is.na(genes)]
        genes = genes[(genes %in% rownames(TPM))]
        exp <- TPM[genes,colnames(temp)]
        exp %<>% rownames_to_column(var = "Gene")
        exp_list[[i]] = left_join(row_tmp,exp, by = "Gene")
        Progress(i, length(sheets))
}
names(exp_list) = sheets
openxlsx::write.xlsx(exp_list, file =  "Yang/GTEx/Integrated-signatures-categorized - 577 GTEx samples ordered.xlsx",
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))


# replace name
