invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#==============
csv_names = paste0("SCT_snn_res.",c(0.01, 0.1, 0.2, 0.5, 0.8, 0.9, 1, 2, 3,4,5))
csv_index = list.files("output/20220418",pattern = ".csv") %>% gsub("_.*","",.) %>% as.integer()
table(1:1661 %in% csv_index)
csv_names = list.files("output/20220418",pattern = ".csv")
deg_list <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(paste0("output/20220418/",csv),row.names = 1)
        tmp = tmp[tmp$p_val_adj < 0.05,]
        tmp$gene = rownames(tmp)
        tmp %<>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
        tmp$resolution = sub("^.*_SCT","SCT",csv) %>% 
                sub(".csv","",.) %>% 
                sub("SCT_snn_res","SCT-snn-res",.) %>%
                sub("_.*","",.) %>%
                sub("SCT-snn-res","SCT_snn_res",.)
        tmp
})

deg = bind_rows(deg_list)
deg %<>% filter(p_val_adj < 0.05)
deg_list = split(deg, f = deg$resolution)
write.xlsx(deg_list, file = paste0(path,"Lung_63_DEG.xlsx"),
           colNames = TRUE, borders = "surrounding")

#===========================================================
csv_names = list.files("output/20220327/Cell_subtype",pattern = ".csv")
deg_list <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(paste0("output/20220327/Cell_subtype/",csv),row.names = 1)
        #tmp = tmp[tmp$p_val_adj < 0.05,]
        tmp %<>% arrange(desc(avg_log2FC))
        tmp$gene = rownames(tmp)
        tmp
})

deg = bind_rows(deg_list)
deg1 = filter(deg, avg_log2FC > 0) %>% group_by(cluster) %>% top_n(50, avg_log2FC)
deg2 = filter(deg) %>% group_by(cluster) %>% top_n(100, abs(avg_log2FC))
deg_list = list(deg1,deg2)
names(deg_list) = c("postive","positve_negative")
write.xlsx(deg_list, file = "output/20220327/WC_30_T_Adj_Dex_vs_Count.xlsx",
           colNames = TRUE, borders = "surrounding")
