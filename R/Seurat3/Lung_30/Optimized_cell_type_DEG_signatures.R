########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("dplyr","cowplot","magrittr","data.table","ggplot2","tidyr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

df_optimization <- readxl::read_excel("doc/20201103_Optimized cell type signatures.xlsx")
df_optimization = df_optimization[!is.na(df_optimization$`Cell type`),]
cellTypes = unique(df_optimization$`Cell type`)

read.path ="Yang/Lung_30/DE_analysis/"

# combine C1+C2+C3 
C_Cell_types_file_name = list.files(path = paste0(read.path,"C_Cell_types/"),
                                    pattern = paste(paste0("_FC0.1_",c("C1","C2","C3"),".csv"),
                                                    collapse = "|"),
                                    full.names = T)
deg_list <- lapply(C_Cell_types_file_name, read.csv, row.names = 1)
C_deg <- bind_rows(deg_list)
C_deg$cluster = "C1.C2.C3"
write.csv(C_deg, file = paste0(read.path,"C_Cell_types/Lung_30-01_FC0.1_C1.C2.C3.csv"))


optimzed_degs_list = list()
for(i in 1:length(cellTypes)){
    # Step 1: Include DEGs of a given cell type logFC >=0.5
    cellType = cellTypes[i]

    C_Cell_types_file_name = list.files(path = paste0(read.path,"C_Cell_types/"),
                                        pattern = paste0("_FC0.1_",cellType,".csv"),
                                        full.names = T)
    degs1 = read.csv(file = C_Cell_types_file_name,row.names = 1)
    markers = degs1[degs1$avg_logFC >=0.5, "gene"]
    # Step 2: Remove the following genes (B_Cell group comparisons)
    remove_groups = which(df_optimization$`Cell type` %in% cellType)
    rm_markers <- c()
    for(k in 1:length(remove_groups)){
        comparison_idx = sub("B","", df_optimization$`Step 2: Remove the following genes(B_Cell group comparisons)`[remove_groups[k]])
        comparison_idx = as.integer(comparison_idx)
        if(comparison_idx < 10) comparison_idx %<>% paste0("0",.)
        B_cell_groups_file_name <- list.files(path = paste0(read.path,"B_Cell_groups"),
                                              pattern = paste0("Lung_30_B_",comparison_idx),
                                              full.names = T)
        #if(length(B_cell_groups_file_name) ==0) next
        degs = read.csv(file = B_cell_groups_file_name,row.names = 1)
        cluster = df_optimization$...5[remove_groups[k]]
        degs = degs[(degs$cluster %in% cluster) & 
                        (degs$avg_logFC > 0) &
                        (degs$p_val_adj < 0.05),]

        rm_markers = c(rm_markers,degs$gene)
        print(paste0('i=',i,"   k=",k))
        print(length(rm_markers))
    }
    purified_markers = markers[!(markers %in% rm_markers)]
    optimzed_degs_list[[i]] = degs1[degs1$gene %in% purified_markers,]
    svMisc::progress(i,max.value = length(cellTypes))
}
optimzed_degs = bind_rows(optimzed_degs_list)
optimzed_degs = optimzed_degs[order(optimzed_degs$cluster, optimzed_degs$p_val_adj),]
rownames(optimzed_degs) = make.unique(optimzed_degs$gene)
write.csv(optimzed_degs, file = paste0(read.path,"Optimized_cell_type_DEG.csv"))
