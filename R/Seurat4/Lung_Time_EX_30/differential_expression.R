invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# if(step == "resolutions"){
res = c(0.5, 1, 2, 5)
csv_names = paste0("SCT_snn_res.",res)


deg_list <- pbapply::pblapply(csv_names, function(csv_name){
        csv = list.files("output/20220524/EX30",pattern = csv_name,full.names = T)
        tmp <- lapply(csv, function(x){
                temp = read.csv(x,row.names = 1)
                temp$gene = rownames(temp)
                temp
        })
        tmp %<>% bind_rows()
        tmp = tmp[tmp$p_val_adj < 0.05,]
        
        tmp %<>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
        tmp
})

names(deg_list) =csv_names

write.xlsx(deg_list, file = paste0(path,"Lung_time_EX30.xlsx"),
           colNames = TRUE, borders = "surrounding")

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



# if(step == "Analysis 2"){# 16~32GB
csv_names = list.files("output/20211022/ASE",pattern = "csv",full.names = T)
DEG <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(csv,row.names = 1) %>% arrange(desc(avg_log2FC))
        tmp$gene = rownames(tmp)
        type = sub(".*[0-9+]-","",csv) %>% sub("_vs_.*","",.)
        tmp$cluster = type
        tmp
}) %>% bind_rows()
write.xlsx(DEG, file = paste0(path,"Lung_SAE_DEG.xlsx"),gridLines = TRUE,
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))



xlsx_list <- list.files(path = "output/20211028",
                        pattern = "2021-10-28-.*xlsx",full.names = T)
deg_list <- pbapply::pblapply(xlsx_list, function(xlsx){
        tmp <- readxl::read_excel(xlsx)
        tmp = tmp[order(tmp$avg_log2FC,decreasing = T), ]
        tmp %<>% filter(p_val_adj < 0.05)
        tmp = tmp[abs(tmp$avg_log2FC) > 0.25, ]
        tmp
})
xlsx_list %<>% gsub(".*2021-10-28-","",.) %>% gsub("\\.xlsx","",.)
names(deg_list) = xlsx_list
write.xlsx(deg_list, file = paste0(path,"Lung_d28_DEG.xlsx"),gridLines = TRUE,
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
