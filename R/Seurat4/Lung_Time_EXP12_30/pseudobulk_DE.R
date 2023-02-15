invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx","Seurat","Libra","MAST","future","Libra",
                   "gplots"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Libra.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#====================================
#==========================
object <- readRDS(file = "data/Lung_time15_SoupX_20230129.rds")
meta.data <- readRDS("output/Lung_time15_metadata_20220523_v2.rds")
meta.data$sample <- meta.data$patient
meta.data$day %<>% factor(levels = paste0("D",c(0,3,7,14,21,28,56,112,122)))
meta.data$orig.ident %<>% factor(levels = c(paste0("EXP12_D",c(0,3,7,14,21,28,56,122)),paste0("EX30_D",c(0,7,14,21,28,56,112))))
if(any(duplicated(colnames(meta.data)))) {
    print(table(duplicated(colnames(meta.data))))
    meta.data[,duplicated(colnames(meta.data))] = NULL
}
if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}
object$All_celles <- "All_celles"
DE <- run_de.1(object, de_family = 'pseudobulk', cell_type_col = "All_celles", 
            label_col = "cell state",replicate_col = "orig.ident",
            ident.1 = "BC-0",min_cells = 3,min_reps = 1,
            de_method = 'edgeR', de_type = 'LRT')
DE_filter <- DE %>% filter(avg_logFC > 0 & p_val_adj < 0.05) %>% 
    group_by(cluster) %>%
    arrange(desc(avg_logFC),.by_group = TRUE)
DE_filter$cell_type = NULL
write.csv(DE_filter,paste0(path,"pseudobulk_de.csv"),row.names = FALSE)

scDE <- readxl::read_excel("output/20230131/Lung_time_EXP12_30_DEGs.xlsx", sheet = "cell state")
scDE %>% filter(cluster == "BC-0") %>% .[["gene"]] %>% head(200) -> scDE_BC0
DE_filter %>% filter(cluster == "BC-0") %>% .[["gene"]] %>% head(200) -> psudoDE_BC0
    
euler_df <- eulerr::euler(list("singlecell_de" = scDE_BC0,
                                   "pseudobulk_de" = psudoDE_BC0))

g <- plot(euler_df, quantities = TRUE, lty = 1:6,  legend = TRUE)
jpeg(paste0(path,"BC0_overlap.jpg"), units="in", width=5, height=3,res=300)
print(g)
dev.off()
