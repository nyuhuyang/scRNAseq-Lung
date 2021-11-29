########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","eulerr","openxlsx",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# =======================================
#if(step == "Analysis 1b"){# 16~32GB
types = c("Cell_subtype","Cell_type")
invitro_DEG <- pbapply::pblapply(types, function(type){
        csv_names = list.files("output/20211019/Analysis 1b",pattern = type,full.names = T)
        csv_names = grep("D28", csv_names, value = T)
        deg <- pbapply::pblapply(csv_names, function(csv, type){
                tmp <- read.csv(csv,row.names = 1) %>% arrange(desc(avg_log2FC))
                tmp$gene = rownames(tmp)
                tmp
        }) %>% bind_rows()
})%>% bind_rows()

rownames(invitro_DEG) =NULL
table(invitro_DEG$p_val_adj <0.05)
table(abs(invitro_DEG$avg_log2FC) >1 )
table(abs(invitro_DEG$avg_log2FC) >0.5 )

#invitro_DEG = invitro_DEG[invitro_DEG$p_val_adj <0.05,]
invitro_DEG = invitro_DEG[abs(invitro_DEG$avg_log2FC) >0.5,]
write.xlsx(invitro_DEG, file = paste0(path,"Lung_invitro_d28_DEG.xlsx"),gridLines = TRUE,
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
# if(step == "Analysis 2"){# 16~32GB
csv_names = list.files("output/20211022/ASE",pattern = "csv",full.names = T)
invivo_DEG <- pbapply::pblapply(csv_names, function(csv){
                                tmp <- read.csv(csv,row.names = 1) %>% arrange(desc(avg_log2FC))
                                tmp$gene = rownames(tmp)
                                type = sub(".*[0-9+]-","",csv) %>% sub("_vs_.*","",.)
                                tmp$type = type
                                tmp
}) %>% bind_rows()

rownames(invivo_DEG) =NULL
table(invivo_DEG$p_val_adj <0.05)
table(abs(invivo_DEG$avg_log2FC) >1 )
table(abs(invivo_DEG$avg_log2FC) >0.5 )

#invivo_DEG = invivo_DEG[invivo_DEG$p_val_adj <0.05,]
invivo_DEG = invivo_DEG[abs(invivo_DEG$avg_log2FC) >0.5,]
write.xlsx(invivo_DEG, file = paste0(path,"Lung_invivo_SAE_DEG.xlsx"),gridLines = TRUE,
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))


#============================
for(fc_value in c(1,0.5,0.25)){
        eulerr(invivo_DEG, key = unique(invivo_DEG$type),
                            group.by = "type", shape =  "circle",#key = c("C1","C2","C3","C4","B_cells"),
                            cut_off = "avg_log2FC", cut_off_value = fc_value,do.print = T,
                            save.path = path, file.name =paste0("Venn_invivo_log2fc_",fc_value,".jpeg"))
}

invivo_DEG %<>% filter(type %in% c("BC","H","IC","p-C","S-Muc","S1","TASC" ))
invivo_DEG$type %<>% paste0("_invivo") 

invitro_DEG$type = gsub("_D28","",invitro_DEG$cluster)
invitro_DEG %<>% filter(type %in% c("BC-m","BC","IC-sq","IC1","S1","TASC"))
invitro_DEG$type %<>% paste0("_invitro") 
invitro_DEG$day_pairs = NULL
invitro_DEG$category = NULL
invitro_DEG = invitro_DEG[,colnames(invivo_DEG)]
DEG = rbind(invivo_DEG, invitro_DEG)
for(fc_value in c(0.5)){
        pos_genes = eulerr(DEG, #key = c("BC","H","IC","p-C","S-Muc","S1","TASC" ),
                            group.by = "type", shape =  "circle",#key = c("C1","C2","C3","C4","B_cells"),
                            cut_off = "avg_log2FC", cut_off_value = fc_value,do.print = T,return.raw = T,
                            save.path = path, file.name =paste0("Venn_vivo_vs_vitro_log2fc_",fc_value,".jpeg"))
}

cell_types = c("TASC","BC","IC", "S1")
Pos_genes_list = list()
for(i in seq_along(cell_types)){
        cell = cell_types[i]
        pos_genes = eulerr(DEG, key = grep(cell,unique(DEG$type),value = T),
               group.by = "type", shape =  "circle",
               cut_off = "avg_log2FC", cut_off_value = 0.5,do.print = F,return.raw = T,
               save.path = path, file.name =paste0("Venn_vivo_vs_vitro_",cell,"_log2fc_1.jpeg"))
        euler_df <- eulerr::euler(pos_genes,shape = "circle")
        pos_genes_list <- as.list(euler_df$original.values)
        Colnames  = paste(names(pos_genes_list), ":",pos_genes_list)
        id <- eulerr:::bit_indexr(8)
        for (k in nrow(id):1) {
                pos_genes_list[[k]] = Reduce(intersect, pos_genes[id[k,]])  %>%
                        setdiff(Reduce(union, pos_genes[!id[k,]]))
        }
        pos_genes_df <- list2df(pos_genes_list)
        colnames(pos_genes_df) = Colnames
        pos_genes_df = pos_genes_df[,sapply(pos_genes_list,length) != 0]
        Pos_genes_list[[i]] = pos_genes_df
        write.xlsx(pos_genes_df, file = paste0(path,"Lung_",cell,"_intersect_DEGs.xlsx"),gridLines = TRUE,
                   colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
        print(cell)
}
names(Pos_genes_list) = cell_types
write.xlsx(Pos_genes_list, file = paste0(path,"Lung_invivo_vitro_SAE_intersect_DEGs.xlsx"),gridLines = TRUE,
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))




