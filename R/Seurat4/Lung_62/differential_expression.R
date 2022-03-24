invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#==============
csv_names = paste0("SCT_snn_res.",c(0.01, 0.1, 0.2, 0.5, 0.8))
csv_index = list.files("output/20220323",pattern = ".csv") %>% gsub("_.*","",.) %>% as.integer()
table(1:376 %in% csv_index)
csv_names = list.files("output/20220323",pattern = ".csv")
deg_list <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(paste0("output/20220323/",csv),row.names = 1)
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
write.xlsx(deg_list, file = paste0(path,"Lung_62_DEG.xlsx"),
           colNames = TRUE, borders = "surrounding")

