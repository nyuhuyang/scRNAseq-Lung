invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# if(step == "resolutions"){
res = c(0.01, 0.05, 0.09, 0.1, 0.5, 0.6, 0.7)
csv_names = paste0("SCT_snn_res.",res)


deg_list <- pbapply::pblapply(csv_names, function(csv){
        csv = list.files("output/20210909",pattern = csv,full.names = T)
        tmp <- read.csv(csv,row.names = 1)
        tmp = tmp[tmp$p_val_adj < 0.05,]
        
        tmp %<>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
        tmp
})

names(deg_list) =csv_names

write.xlsx(deg_list, file = paste0(path,"Lung_time6_DEG.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#if(step == "Analysis 1a"){# 16~32GB
types = c("Cell_subtype","Cell_type")
DEG1 <- pbapply::pblapply(types, function(type){
        csv_names = list.files("output/20211018/Analysis 1a",pattern = type,full.names = T)
        deg <- pbapply::pblapply(csv_names, function(csv, type){
                tmp <- read.csv(csv,row.names = 1) %>% arrange(desc(avg_log2FC))
                tmp$gene = rownames(tmp)
                tmp
        }) %>% bind_rows()
})%>% bind_rows()

rownames(DEG1) =NULL
#write.xlsx(DEG1, file = paste0(path,"Lung_time6_DEG_Analysis_1a.xlsx"),
#           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))


#if(step == "Analysis 1b"){# 16~32GB
types = c("Cell_subtype","Cell_type")
DEG2 <- pbapply::pblapply(types, function(type){
        csv_names = list.files("output/20211019/Analysis 1b",pattern = type,full.names = T)
        deg <- pbapply::pblapply(csv_names, function(csv, type){
                tmp <- read.csv(csv,row.names = 1) %>% arrange(desc(avg_log2FC))
                tmp$gene = rownames(tmp)
                tmp
        }) %>% bind_rows()
})%>% bind_rows()

rownames(DEG2) =NULL
DEG_list = list(DEG1,DEG2)
names(DEG_list) = c("Analysis 1a","Analysis 1b")
write.xlsx(DEG_list, file = paste0(path,"Lung_time6_DEG_Analysis_1.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))


#if(step == "Analysis 2") {# 16~32GB
types = c("D0", "D3","D7","D14","D21", "D28")
DEG <- pbapply::pblapply(types, function(type){
        csv_names = list.files("output/20211021/Analysis 2",pattern = type,full.names = T)
        deg <- pbapply::pblapply(csv_names, function(csv, type){
                tmp <- read.csv(csv,row.names = 1) %>% arrange(desc(avg_log2FC))
                tmp$gene = rownames(tmp)
                tmp
        }) %>% bind_rows()
})%>% bind_rows()

rownames(DEG) =NULL
table(DEG$p_val_adj <0.05)
table(abs(DEG$avg_log2FC) >1 )
DEG = DEG[DEG$p_val_adj <0.05,]
DEG = DEG[abs(DEG$avg_log2FC) >1,]

DEG_list = split(DEG, f = DEG$type)
write.xlsx(DEG_list, file = paste0(path,"Lung_time6_DEG_Analysis_2.xlsx~"),gridLines = TRUE,
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
