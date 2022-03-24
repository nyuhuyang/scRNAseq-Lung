invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots","kableExtra"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

 ############### step == "resolutions ############### 

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

############### step == "cell_types" ############### 
opts = data.frame(ident = c(rep("Cell_subtype",49),
                            rep("Cell_type",31),
                            rep("Cell_group",26), #81~106
                            rep("UMAP_land",20),
                            rep("Family",7),
                            rep("Superfamily",3)),
                  num = c(1:49,
                          1:31,
                          1:26,
                          1:20,
                          1:7,
                          1:3)
)
Cell_category = unique(opts$ident)
deg_list <- list()
for(i in seq_along(Cell_category)){
        csv_names = list.files("output/20220216/cell_types",pattern = Cell_category[i],full.names = T)
        all_idx = which(opts$ident %in% Cell_category[i])
        idx <- gsub("output/20220216/cell_types/","",csv_names) %>% gsub("-.*","",.) %>% as.integer()
        print(table(all_idx %in% idx))
        all_idx[!all_idx %in% idx]
        print(paste(Cell_category[i], "missing",all_idx[!(all_idx %in% idx)]))
        deg <- pbapply::pblapply(csv_names, function(csv){
                tmp <- read.csv(csv,row.names = 1)
                tmp$gene = rownames(tmp)
                tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
                if(tmp$cluster == TRUE) tmp$cluster = "T"
                return(tmp)
        }) %>% bind_rows
        deg = deg[deg$p_val_adj < 0.05,]
        deg_list[[i]] = deg
}
names(deg_list) =Cell_category

write.xlsx(deg_list, file = paste0(path,"Lung_30_DEG_Cell.category.xlsx"),
           colNames = TRUE, borders = "surrounding")

# On Sep 24, 2021, at 13:28, Renat Shaykhiev <res2003@med.cornell.edu> wrote:
# DEG analysis 1 - between TASC and related cell populations (with volcanos; cells in the same group from different samples mixed):
save.path = "Yang/Lung_30/hg38/DE_analysis/"
csv_names <- list.files("output/20211230/DEG analysis 1")
deg <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(paste0("output/20211230/DEG analysis 1/", csv),row.names = 1)
        tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
        region = sub(".*_","",csv) %>% sub(".csv","",.)
        tmp$gene = rownames(tmp)
        tmp$region = switch(region,
                            "D" = "distal",
                            "DT" = "distal+terminal")
        return(tmp)
}) %>% bind_rows() %>% filter(p_val_adj < 0.05)

deg_list = split(deg, f = deg$region)
write.xlsx(deg_list, file = paste0(save.path,"Lung_30_DEG_TASC_related.xlsx"),
           colNames = TRUE, borders = "surrounding")

# DEG analysis 2 – between different subsets of TASCs
csv_names <- list.files("output/20211230/DEG analysis 2")
deg <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(paste0("output/20211230/DEG analysis 2/", csv),row.names = 1)
        tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
        return(tmp)
}) %>% bind_rows() %>% filter(p_val_adj < 0.05)


write.xlsx(deg, file = "output/20211230/DEG analysis 2/Lung_30_DEG_TASC_subset.xlsx",
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#DEG analysis 3 - between sample groups
#Option 1 – considering all cells in the group together
Cell_types <- c("Cell_subtype","Cell_type","UMAP_land","Family","Superfamily")
meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")

categories = lapply(Cell_types, function(Cell_type){
        unique(meta.data[,Cell_type])
})
names(categories) = Cell_types
all_cells = unique(unlist(categories,use.names = F))
all_cells = all_cells[all_cells != "Un"]

csv_names = list.files("output/20211012/DEG analysis 3-option 1",pattern = ".csv",full.names = T)
idx <- gsub("output/20211012/DEG analysis 3-option 1/","",csv_names) %>% gsub("-.*","",.) %>% as.integer()
all_idx = 1:900
table(all_idx %in% idx)
all_idx[!all_idx %in% idx]
print(paste("missing",sort(as.character(all_idx[!(all_idx %in% idx)]))))
missing_idx = paste(all_idx[!all_idx %in% idx],collapse = ",")
write.csv(missing_idx,file = paste0(path,"missing_idx.csv"))

deg <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(csv,row.names = 1)
        tmp$gene = rownames(tmp)
        tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
        if(tmp$type == TRUE) tmp$type = "T"
        return(tmp)
}) %>% bind_rows
deg = deg[deg$p_val_adj < 0.05,]
rownames(deg) = NULL
deg_list = split(deg,f = deg$type)

write.xlsx(deg_list, file = paste0(path,"Lung_30_DEG_between_groups_option1.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))


#DEG analysis 3 - between sample groups
#Option 2 – considering using average expression
Cell_types <- c("Cell_subtype","Cell_type","UMAP_land","Family","Superfamily")
meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")

categories = lapply(Cell_types, function(Cell_type){
        unique(meta.data[,Cell_type])
})
names(categories) = Cell_types
all_cells = unique(unlist(categories,use.names = F))
all_cells = all_cells[all_cells != "Un"]

csv_names = list.files("output/20211012/DEG analysis 3-option 2",pattern = ".csv",full.names = T)
idx <- gsub("output/20211012/DEG analysis 3-option 2/","",csv_names) %>% gsub("-.*","",.) %>% as.integer()
all_idx = 1:900
table(all_idx %in% idx)
all_idx[!all_idx %in% idx]
missing_idx = paste(all_idx[!all_idx %in% idx],collapse = ",")
write.csv(missing_idx,file = paste0(path,"missing_idx.csv"))
print(paste("missing",sort(as.character(all_idx[!(all_idx %in% idx)]))))

deg <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(csv,row.names = 1)
        tmp$gene = rownames(tmp)
        tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
        if(tmp$type == TRUE) tmp$type = "T"
        return(tmp)
}) %>% bind_rows
deg = deg[deg$p_val < 0.05,]
rownames(deg) = NULL
deg_list = split(deg,f = deg$type)

write.xlsx(deg_list, file = paste0(path,"Lung_30_DEG_between_groups_option2.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

csv_names = csv_names[grep("COPD",csv_names)]
deg <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(csv,row.names = 1)
        tmp$gene = rownames(tmp)
        tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
        if(tmp$type == TRUE) tmp$type = "T"
        return(tmp)
}) %>% bind_rows
deg = deg[deg$p_val < 0.05,]
table(deg$type,deg$DE_pairs) %>% kable %>% kable_styling()
deg = deg[deg$avg_log2FC > 0.5,]

deg_list = split(deg,f = deg$type)

#=========step == "TASCs"============================
csv_names <- list.files("output/20211129",pattern = "_TACS.csv",full.names = T)
deg <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(csv,row.names = 1)
        tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
        return(tmp)
}) %>% bind_rows
deg = deg[deg$p_val < 0.05,]
rownames(deg) = NULL
deg_list = split(deg,f = deg$cluster)

write.xlsx(deg_list, file = paste0(path,"TACS_DEG_SCGB1A1_SCGB3A2.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))



csv_names <- list.files("output/20211129",pattern = "SCGB1A1 high vs low at.*_TACS.csv",full.names = T)
deg <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(csv,row.names = 1)
        tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
        return(tmp)
}) %>% bind_rows
deg = deg[deg$p_val < 0.05,]
deg[grep("SFTPA2",deg$gene),c("avg_log2FC","cluster")]
deg[grep("SFTPA1",deg$gene),c("avg_log2FC","cluster")]
deg[grep("HOPX",deg$gene),c("avg_log2FC","cluster")]

rownames(deg) = NULL
deg_list = split(deg,f = deg$cluster)

write.xlsx(deg_list, file = paste0(path,"TACS_DEG_SCGB1A1.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#=========step == "ASE"============================
save.path = "Yang/Lung_30/hg38/DE_analysis/"
csv_names <- list.files("output/20220210/ASE",full.names = T)
deg <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(csv,row.names = 1)
        tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
        cell_type = sub("_vs_.*","",csv) %>% sub(".*[0-9+]-","",.)
        region = sub(".*_","",csv) %>% sub(".csv","",.)
        tmp$gene = rownames(tmp)
        tmp$cell_type = cell_type
        tmp$region = switch(region,
                            "D" = "distal",
                            "DT" = "distal+terminal")
        return(tmp)
}) %>% bind_rows() %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0)

deg_list = split(deg, f = deg$region)
write.xlsx(deg_list, file = paste0(save.path,"Lung_30_DEG_ASE.xlsx"),
           colNames = TRUE, borders = "surrounding",overwrite = T)

#========="deconvolution"============================
#Level 3 ASE: BC+IC, S (S1+S-Muc+TASC), C (C1+C-s), p-C, Ion, NE
#Level 3 ASE+AT: all ASE above + AT
#Level 3 Structural: Cr, Fb (Fb1+2+3+4), SM+Pr (SM1+2+3+Pr), En-v, En-c (En-c1+En-ca)
#Level 3 Immune: B, PC, T, MC, Neu, MPS (Mon+M1+M1-2+M2+c-DC)  
#Level 3 combined: all above combined
names_list <- list("Level 3 ASE" = c("BC-IC","S","C","Ion","NE"),
                   "Level 3 ASE+AT" = c("BC-IC","S","C","Ion","NE","AT"),
                   "Level 3 Structural" = c("Cr","Fb","SM+Pr","En-v","En-c"),
                   "Level 3 Immune" = c("B", "PC", "T", "MC", "Neu", "MPS"))
names_list[["Level 3 combined"]] = unlist(names_list,use.names = FALSE)
save.path = "Yang/Lung_30/hg38/DE_analysis/"
categories = c("Cell_subtype","Cell_type","Cell_group","UMAP_land","Family","Superfamily")

deg_list  = pbapply::pblapply(categories, function(category) {
        readxl::read_excel(paste0("output/20220216/","Lung_30_DEG_Cell.category_FC1.xlsx"),sheet = category)
})
names(deg_list) = categories
deg = deg_list[["Cell_group"]]
deg1 = deg_list[["Cell_subtype"]]
deg %<>% rbind(deg1[deg1$cluster %in% "p-C",])
deg_list1  = pbapply::pblapply(names_list, function(clusters) {
        deg %>% filter(cluster %in% clusters) %>%  filter(abs(avg_log2FC) > 1)
        })
deg_list %<>% c(deg_list1)
openxlsx::write.xlsx(deg_list, file =  paste0("output/20220216/","Lung_30_DEG_Cell.category_neg_posLog2FC1.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",overwrite = T)

deg_list  = pbapply::pblapply(deg_list, function(deg) {
        deg %>% filter(avg_log2FC > 1)
})
openxlsx::write.xlsx(deg_list, file =  paste0(path,"Lung_30_DEG_Cell.category_posLog2FC1.xlsx"),
                     colNames = TRUE,row.names = F,borders = "surrounding",overwrite = T)
