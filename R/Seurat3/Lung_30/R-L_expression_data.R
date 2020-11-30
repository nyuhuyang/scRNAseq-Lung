# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","magrittr",
                   "openxlsx","tidyr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#======1.1 read template =========================
R_L_template <- readxl::read_excel("doc/R-L expression data per cell type for each sample.xlsx",
                                   sheet = "each sample")
gene_R = R_L_template$R
gene_L = R_L_template$L
genes = unique(c(gene_R,gene_L))
colnames(R_L_template) %<>% gsub("\\..*","",.)
R_L_template = as.data.frame(R_L_template[,unique(colnames(R_L_template))])
cell.types = unique(colnames(R_L_template[,-c(1:3)]))
genes[!(genes %in% rownames(object[["RNA"]]@data))]

changeName <- data.frame("Aliases" = c("ACVR1A","Il13RA1", "IL36R", "TNFRSF3", "TSLPR", "BMP9", "VEGFD"),
                         "gene" = c("ACVR1", "IL13RA1", "IL36RN", "LTBR", "CRLF2", "GDF2", "FIGF"))
for(i in 1:nrow(changeName)){
    gene_R %<>% gsub(changeName$Aliases[i], changeName$gene[i], .)
    gene_L %<>% gsub(changeName$Aliases[i], changeName$gene[i], .)
}
genes = unique(c(gene_R,gene_L))
genes[!(genes %in% rownames(object[["RNA"]]@data))]

#======1.2 load  Seurat =========================
# load files
object = readRDS(file = "data/Lung_30_20200710.rds") 
DefaultAssay(object) = "RNA"
Idents(object) = "annotations3"
object %<>% subset(idents = cell.types)
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object)= "orig.ident"

samples = levels(Idents(object))
R_L_exp_list <- vector(mode = "list", length = length(samples))
for(j in 1:length(samples)){
    sub_object = subset(object, idents = samples[j])
    Idents(sub_object) = "annotations3"
    sub_object %<>% sortIdent()
    R_L_exp = AverageExpression(sub_object,assays = "RNA",verbose = T,
                                features = unique(c(gene_R,gene_L)))
    tempR <- R_L_template
    tempL <- R_L_template[,-c(1:3)]
    tempR[,colnames(R_L_exp$RNA)] = R_L_exp$RNA[gene_R,]
    tempL[,colnames(R_L_exp$RNA)] = R_L_exp$RNA[gene_L,]
    temp <- cbind(tempR,tempL)
    temp[is.na(temp)] = 0
    R_L_exp_list[[j]]=temp
    svMisc::progress(j, max.value = length(samples))
}
names(R_L_exp_list) = samples
write.xlsx(R_L_exp_list, file = paste0(path,"Lung_30_R-L expression data per cell type for each sample.xlsx"),
           colNames = TRUE, rowNames = FALSE,borders = "surrounding",
           colWidths = c(NA, "auto", "auto"),
           firstActiveCol = 4)

# Prepare final Output R-L cell-cell interactions per sample
samples <- c("UNC-48-P","UNC-55-P","UNC-66-P","UNC-69-P","UNC-71-P","UNC-48-D",
             "UNC-55-D","UNC-66-D","UNC-69-D","UNC-71-D","UNC-44-D","UNC-54-D",
             "UNC-57-D","UNC-67-D","UNC-70-D","CU-12-D","CU-12-D-R","VU-27-D",
             "UNC-48-T","UNC-55-T","UNC-66-T","UNC-69-T","UNC-71-T","UNC-44-T",
             "CU-12-T","UNC-51-D","UNC-52-D","UNC-61-D","VU-19-D","VU-37-D")
df_list <- list()
for(i in 1:length(samples)) {
    s = samples[i]
    temp <- readxl::read_excel("Yang/Lung_30/Cell_Phone_DB/R-L interaction score-samples.xlsx",
                                       sheet = s)
    colnames(temp)[1] = "index"
    temp %<>% gather(key = "Cell-cell pair", value = "value", 
                     -c("index","R","L","R-L pair"))
    colnames(temp)[grep("value",colnames(temp))]= s
    temp = temp[order(temp$index),]
    df_list[[i]]= temp
    svMisc::progress(i, length(samples))
}

sapply(df_list, function(x) identical(temp[,1:5],x[,1:5]))

df <- temp[,1:5]
for(i in 1:length(samples)) {
    s = samples[i]
    df %<>% cbind(df_list[[i]][,s])
    svMisc::progress(i, length(samples))
}
colnames(df)[1] = ""
write.xlsx(df, file = "Yang/Lung_30/Cell_Phone_DB/R-L interaction score-summary.xlsx",
           colNames = TRUE, rowNames = FALSE,borders = "surrounding",
           colWidths = c(NA, "auto", "auto"),
           firstActiveCol = 6)
# ===========Statistics output =============
group_list <- list("P" = c("UNC-48-P","UNC-55-P","UNC-66-P","UNC-69-P","UNC-71-P"),
                   "D" = c("UNC-48-D","UNC-55-D","UNC-66-D","UNC-69-D","UNC-71-D","UNC-44-D","UNC-54-D",
                           "UNC-57-D","UNC-67-D","UNC-70-D","CU-12-D","CU-12-D-R","VU-27-D"),
                   "P+D" = c("UNC-48-P","UNC-55-P","UNC-66-P","UNC-69-P","UNC-71-P",
                             "UNC-48-D","UNC-55-D","UNC-66-D","UNC-69-D","UNC-71-D","UNC-44-D","UNC-54-D",
                             "UNC-57-D","UNC-67-D","UNC-70-D","CU-12-D","CU-12-D-R","VU-27-D"),
                   "T" = c("UNC-48-T","UNC-55-T","UNC-66-T","UNC-69-T","UNC-71-T","UNC-44-T",
                           "CU-12-T"),
                   "D+T" = c("UNC-48-D","UNC-55-D","UNC-66-D","UNC-69-D","UNC-71-D","UNC-44-D","UNC-54-D",
                             "UNC-57-D","UNC-67-D","UNC-70-D","CU-12-D","CU-12-D-R","VU-27-D",
                             "UNC-48-T","UNC-55-T","UNC-66-T","UNC-69-T","UNC-71-T","UNC-44-T",
                           "CU-12-T"),
                   "P+T" = c("UNC-48-P","UNC-55-P","UNC-66-P","UNC-69-P","UNC-71-P",
                           "UNC-48-T","UNC-55-T","UNC-66-T","UNC-69-T","UNC-71-T","UNC-44-T",
                           "CU-12-T"),
                   "COPD" = c("UNC-51-D","UNC-52-D","UNC-61-D","VU-19-D","VU-37-D"))
comp <- c("P vs D","P vs T","P vs D+T","D vs P","D vs T","D vs P+T","D vs COPD",
          "T vs P","T vs D","T vs P+D","COPD vs D")
df <- readxl::read_excel(path = "Yang/Lung_30/Cell_Phone_DB/R-L interaction score-summary.xlsx")
df_res <- df[,1:5]
res <- sapply(group_list,function(x) rowSums(df[,x]))
df_res %<>% cbind(res)

keep = rowSums(res) >0
df = df[keep,]
df_res = df_res[keep,]

res <- sapply(comp,function(x) {
    comp1 = sub(" vs.*","",x)
    comp2 = sub(".*vs ","",x)
    log2(df_res[,comp1]+1)-log2(df_res[,comp2]+1)
})
df_res %<>% cbind(res)
# wilcox.test(y1,y2,paired=TRUE)
res <- sapply(comp,function(s) {
    comp1 = sub(" vs.*","",s)
    comp2 = sub(".*vs ","",s)
    cells.1 = group_list[[comp1]]
    cells.2 = group_list[[comp2]]
    
    group.info <- data.frame(row.names = c(cells.1, cells.2))
    group.info[cells.1, "group"] <- "Group1"
    group.info[cells.2, "group"] <- "Group2"
    group.info[, "group"] <- factor(x = group.info[, "group"])
    
    data.use <- df[, rownames(x = group.info)]
    p_val <- pbapply::pbsapply(
        X = 1:nrow(x = data.use),
        FUN = function(x) {
            return(wilcox.test(unlist(data.use[x, ]) ~ group.info[, "group"], 
                               alternative = c("two.sided"))$p.value)
        }
    )
    
})

p_res <- res
p_res[is.nan(p_res)] = 1
df_res %<>% cbind(p_res)

p_val_adj <- pbapply::pbapply(X = p_res,
                              MARGIN = 2, 
                              FUN = function(x) p.adjust(
                                   p = x,
                                   method = "fdr",
                                   n = length(x)))
df_res %<>% cbind(p_val_adj)
openxlsx::write.xlsx(df_res, 
                     file = paste0(path,"Statistics for R-L cell-cell interactions per groups.xlsx"),
                     colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
