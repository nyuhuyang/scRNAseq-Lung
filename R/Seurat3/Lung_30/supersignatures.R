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
        read.csv(x,row.names = 1) %>% 
                mutate(pct.1_pct.2 = pct.1-pct.2) %>%
                filter(avg_logFC >=1) %>%
                filter(pct.1 >=0.5) %>%
                filter(pct.2 < 0.3) %>%
                filter(pct.1_pct.2 > 0.4)
        
        
})
supersignatures = bind_rows(deg_list)
rownames(supersignatures) = make.unique(supersignatures$gene)
openxlsx::write.xlsx(supersignatures, 
                     file =  paste0(paste0(read.path,"C_Cell_types/supersignatures.xlsx")),
                     colNames = TRUE,row.names = T,borders = "surrounding",colWidths = c(NA, "auto", "auto"))
