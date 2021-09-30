# conda activate r4.0
library(dplyr)
library(magrittr)
library(stringr)
library(ggsci)
library(enrichR)
library(openxlsx)
library(progress)
library(pbapply)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

set.seed(101)
read.path = "Yang/Lung_30/hg19/DE_analysis/"

# select pathway
wb <- loadWorkbook("Yang/Lung_30/hg38/GSEA/Enrichr/Enrichr libraries.xlsx")
df <- readxl::read_excel("Yang/Lung_30/hg38/GSEA/Enrichr/Enrichr libraries.xlsx")

# filter by color
df_color <- data.frame()
for(i in 1:length(wb$styleObjects)){
    x = wb$styleObjects[[i]]
    if(is.null(x$style$fill)) next
    if(x$style$fill$fillFg == "FFFFFF00"){
        df_color = rbind(df_color, data.frame(ri =x$rows, ci = x$cols, col = "FFFFFF00"))
    }
}
dbs = df[unique(df_color$ri)-1, "Gene-set Library"] %>% pull %>% sort


#====cell type-specific gene ==============

save.path = "Yang/Lung_30/hg38/GSEA/Enrichr/"
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

deg = readxl::read_excel("Yang/Lung_30/hg38/DE_analysis/Lung_30_DEG_Cell.category.xlsx",
                   sheet = "Cell_subtype")
deg[deg$p_val == 0,"p_val"] = min(deg[deg$p_val > 0,"p_val"], .Machine$double.xmin)
deg$log2FC_log10p = deg$avg_log2FC*(-log10(deg$p_val))
cell.type_list <- list("Epi" = c("BC","IC","S1","TASC", "H", "p-C", "C1","C-s", "Ion", "NE","ME",
                                 "G-Ser","G-Muc","S-Muc","AT1","AT2"),
                       "Str" = c("Fb1","Fb2","Fb3","Fb4","Cr","Gli","SM1","SM2","SM3",
                                 "Pr", "En-a","En-c1","En-ca","En-l","En-SM","En-v"),
                       "Im" = c("MC", "Neu", "Mon", "M1", "M1-2","M2","cDC","pDC",
                                "B","PC", "CD4-T1","CD4-T-ifn","CD8-T1", "CD8-T-NK","NK"))

cell.types = unlist(cell.type_list,use.names = F)
deg_list <- split(deg,f = deg$cluster)

enrichedRes <- list()
for(i in 1:length(cell.types)){
    cell = cell.types[i]
    geneRank =  deg_list[[cell.types[i]]]
    colnames(geneRank) %<>% gsub(paste0(cell,"_"),"",.)
    geneRank = geneRank[order(geneRank["log2FC_log10p"]),c("gene","avg_log2FC")]  %>% tibble::deframe()
    tmp <- enrichr(names(geneRank), dbs) #dbs[-which(dbs %in% "dbGaP")])
    # record and remove empty element in tmp
    emp <- c()
    for(k in seq_along(tmp)) {
        if(nrow(tmp[[k]]) > 0 ) {
            tmp[[k]][,"library"] = names(tmp[k])
        } else emp = c(emp, k)
    }
    if(!is.null(emp)) tmp[emp] = NULL
    
    tmp = bind_rows(tmp)
    tmp = tmp[tmp$Adjusted.P.value < 0.05,]
    if(nrow(tmp) > 0 ) {
        tmp$cell.types = cell.types[i]
    } else tmp =NULL
    enrichedRes[[i]] = tmp
    Progress(i, length(cell.types))
}
df_enrichedRes =  bind_rows(enrichedRes)
df_enrichedRes = df_enrichedRes[df_enrichedRes$Adjusted.P.value < 0.05, ]
write.xlsx(df_enrichedRes, asTable = F,
           file = paste0(save.path,"enrichR_DEG_Cell.category.xlsx"),
           borders = "surrounding")
df_enrichedRes <- readxl::read_excel(paste0(save.path,"enrichR_EVG_all.xlsx"))

enrichedRes_list <- split(df_enrichedRes, f = df_enrichedRes$library)

for(i in seq_along(enrichedRes_list)){
    write.csv(enrichedRes_list[[i]], file = paste0(save.path,"enrichR_",
                                                   names(enrichedRes_list)[i],".csv"),
              row.names = FALSE)
    Progress(i, length(enrichedRes_list))
    
}
names(enrichedRes_list) %<>% substr(1, 30)   
write.xlsx(enrichedRes_list, asTable = F,
           file = paste0(save.path,"enrichR_celltypes.xlsx"),
           borders = "surrounding")
#how many samples (and how many different tissues) were used in this HGA data set? 
library(data.table)
HGA <- fread('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Human_Gene_Atlas',
             fill=TRUE,comment.char=.)
# For GTEx, was the v8 release of GTEx included in Enrichr? This will give an idea how many samples were included overall.
GTEx <- fread('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GTEx_Tissue_Sample_Gene_Expression_Profiles_up',
              fill=TRUE)
GTEx <- fread("/Users/yah2014/Downloads/GTEx_Tissue_Sample_Gene_Expression_Profiles_up.txt",fill=TRUE,
              sep="\t")
