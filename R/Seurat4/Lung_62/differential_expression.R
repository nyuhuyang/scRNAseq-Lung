invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#==============
csv_names = paste0("RNA_snn_res.",c(0.01, 0.1, 0.2, 0.5))

deg_list <- pbapply::pblapply(csv_names, function(csv){
        csv = list.files("output/20210912",pattern = csv,full.names = T)
        tmp <- read.csv(csv,row.names = 1)
        tmp = tmp[tmp$p_val_adj < 0.05,]
        
        tmp %<>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
        tmp
})

names(deg_list) = csv_names

write.xlsx(deg_list, file = paste0(path,"Lung_62_DEG.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

