invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))


path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

res = c("UMAP_res=0.8","UMAP_res=1","UMAP_res=2","UMAP_res=3","UMAP_res=4","UMAP_res=5")
opts = data.frame(ident = c(rep("UMAP_res=0.8",84),
                            rep("UMAP_res=1",93),
                            rep("UMAP_res=2",137),
                            rep("UMAP_res=3",174),
                            rep("UMAP_res=4",200),
                            rep("UMAP_res=5",227)),
                  num = c(0:83,
                          0:92,
                          0:136,
                          0:173,
                          0:199,
                          0:226)
)

deg_list <- list()
for(i in seq_along(res)){
        csv_names = list.files("output/20210903",pattern = res[i],full.names = T)
        all_idx = which(opts$ident %in% res[i])
        idx <- gsub("output/20210903/","",csv_names) %>% gsub("_.*","",.) %>% as.integer()
        table(all_idx %in% idx)
        print(paste(res[i], "missing",all_idx[!(all_idx %in% idx)]))
        deg <- pbapply::pblapply(csv_names, function(csv){
                tmp <- read.csv(csv,row.names = 1)
                tmp$gene = rownames(tmp)
                tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
                tmp
        }) %>% bind_rows
        deg = deg[deg$p_val_adj < 0.05,]
        deg_list[[i]] = deg
}
names(deg_list) =res

write.xlsx(deg_list, file = paste0(path,"Lung_30_DEG_20UMAP.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

