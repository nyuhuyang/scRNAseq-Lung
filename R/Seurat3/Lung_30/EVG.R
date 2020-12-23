########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("dplyr","magrittr","data.table","pbapply","tibble","tidyr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#=======
Expression = data.table::fread(file = "output/20201215/Lung_30-Cell.types_expression.txt",
                   sep = "\t") %>% as.data.frame()
Expression %<>% column_to_rownames(var = "V1")
Expression[1:4,1:4]

anno <- readxl::read_excel("doc/Annotations/Cell type abbreviation.xlsx")

df <- readxl::read_excel("doc/Chord diagram cell order - updated abbreviations 12-14-20.xlsx",col_names = T)
superfamily <- c("Epithelial","Structural","Immune")
EVG <- vector(mode = "list",length = length(superfamily))
for(i in seq_along(superfamily)){
    selected.cell.types <- df$cell.types[df$Cell.group %in% superfamily[i]] %>% sort
    gene_list = pbsapply(selected.cell.types,function(ident.1){
        ident.2 = selected.cell.types[!(selected.cell.types %in% ident.1)]
        keep1 = Expression[,paste0(ident.1,".pct")] >= 1/3
        keep2 = Expression[,ident.1]/rowMeans(Expression[,ident.2]) >= 2
       rownames(Expression)[keep1 & keep2]
    })
    gene_list %<>% list2df
    gene_list[is.na(gene_list)] = ""
    EVG[[i]] = gene_list
}
names(EVG) = superfamily
openxlsx::write.xlsx(EVG, file =  paste0(path,"Lung_30-EVGs.xlsx"),
                     colNames = TRUE,row.names = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))


#===== generate FC===================
EVG_df <- pblapply(EVG,function(x) {
    temp = gather(as.data.frame(x),cell.type,genes)
    temp[temp$genes != "",]
    }) %>% bind_rows

res_list <- vector(mode = "list",length = length(superfamily))
cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix("", n-nrow(x), ncol(x))))) 
}

for(i in seq_along(superfamily)){
    selected.cell.types <- df$cell.types[df$Cell.group %in% superfamily[i]] %>% sort
    genes.de = pblapply(selected.cell.types,function(ident.1){
        ident.2 = selected.cell.types[!(selected.cell.types %in% ident.1)]
        genes  = EVG_df %>% filter(cell.type == ident.1) %>% dplyr::select(genes) %>% pull
        
        keep1 = Expression[,paste0(ident.1,".pct")] >= 1/3
        keep2 = Expression[,ident.1]/rowMeans(Expression[,ident.2]) >= 2
        rownames(Expression)[keep1 & keep2]
        de = data.frame(gene  = genes,
                        exp.1 = Expression[genes,ident.1],
                        exp.2 = rowMeans(Expression[genes,ident.2]),
                        pct.1 = Expression[genes,paste0(ident.1,".pct")],
                        pct.2 = rowMeans(Expression[genes,paste0(ident.2,".pct")]),
                        FC = Expression[genes,ident.1]/rowMeans(Expression[genes,ident.2])
        )
        de %<>% arrange(desc(FC))
        colnames(de) %<>% paste0(ident.1,"_",.)
        return(de)
    })
    res <- do.call(cbind.fill, genes.de)
    rownames(res) = NULL
    res_list[[i]] = res
}

names(res_list) = superfamily
openxlsx::write.xlsx(res_list, file =  paste0(path,"Lung_30-EVGs-full.xlsx"),
                     colNames = TRUE,row.names = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

