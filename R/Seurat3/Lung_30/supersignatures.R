invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
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

df_TPM = fread("data/RNA-seq/GTEx-Lung-tpm.csv",header = T)
if(df_TPM$V1[1] == "Age.Bracket") df_TPM = df_TPM[-1,]
TPM <- column_to_rownames(df_TPM, var = "V1") %>% sapply(as.numeric) %>%
        as.data.frame
TPM = cbind("V1"=df_TPM$V1,TPM)
rownames(TPM) = TPM$V1

temp = readxl::read_excel("Yang/GTEx/Cell type EVG genes - 577 GTEx samples ordered.xlsx",sheet = "AT1",)
colnames(temp)[1] = "V1"
exp_list <- list()
for(i in seq_along(cell.types)) {
        genes <- supersignatures[supersignatures$cluster == cell.types[i],"gene"]
        genes = genes[genes != ""]
        genes = genes[!is.na(genes)]
        exp <- TPM[genes,colnames(temp)]
        genes = genes[genes %in% exp$V1]
        exp = exp[match(genes, exp$V1),]
        exp_list[[i]] = exp#rbind.data.frame(temp,exp)
        colnames(exp_list[[i]])[1] = ""
        Progress(i, length(cell.types))
}

names(exp_list) = cell.types
openxlsx::write.xlsx(exp_list, file =  "Yang/GTEx/Cell type supersignatures genes - 577 GTEx samples ordered.xlsx",
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))
