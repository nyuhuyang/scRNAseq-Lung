invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))


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
                            rep("UMAP_land",20),
                            rep("Family",7),
                            rep("Superfamily",3)),
                  num = c(1:49,
                          1:31,
                          1:20,
                          1:7,
                          1:3)
)
Cell_category = unique(opts$ident)
deg_list <- list()
for(i in seq_along(Cell_category)){
        csv_names = list.files("output/20211006/cell_types",pattern = Cell_category[i],full.names = T)
        all_idx = which(opts$ident %in% Cell_category[i])
        idx <- gsub("output/20211006/cell_types/","",csv_names) %>% gsub("-.*","",.) %>% as.integer()
        print(table(all_idx %in% idx))
        #all_idx[!all_idx %in% idx]
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
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

# On Sep 24, 2021, at 13:28, Renat Shaykhiev <res2003@med.cornell.edu> wrote:
# DEG analysis 1 - between TASC and related cell populations (with volcanos; cells in the same group from different samples mixed):
save.path = "/Yang/Lung_30/hg38/DE_analysis/"
xlsx_names <- list.files("output/20210924", pattern = "2021-09-24-")
deg <- pbapply::pblapply(xlsx_names, function(xlsx){
        tmp = readxl::read_excel(paste0("output/20210924/",xlsx))
        tmp[tmp$p_val_adj < 0.05,]
        return(tmp)
})
names(deg) = gsub("2021-09-24-","",xlsx_names) %>% gsub("\\.xlsx$","",.)
write.xlsx(deg, file = paste0(path,"DEG_analysis_1/Lung_30_DEG_TASC_related.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

# DEG analysis 2 – between different subsets of TASCs
xlsx_names <- list.files("output/20210927", pattern = ".xlsx",)
deg <- pbapply::pblapply(xlsx_names, function(xlsx){
        tmp = readxl::read_excel(paste0("output/20210927/",xlsx))
        tmp[tmp$p_val_adj < 0.05,]
        return(tmp)
})
names(deg) = gsub("2021-09-28-","",xlsx_names) %>% gsub("\\.xlsx$","",.)
write.xlsx(deg, file = paste0(path,"Lung_30_DEG_TASC_subset.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#DEG analysis 3 - between sample groups
#Option 1 – considering all cells in the group together
Cell_types <- c("Cell_subtype","Cell_type","Family","Superfamily")
categories = lapply(Cell_types, function(Cell_type){
        unique(object@meta.data[,Cell_type])
})
names(categories) = Cell_types
all_cells = unique(unlist(categories,use.names = F))
all_cells = all_cells[all_cells != "Un"]

csv_names = list.files("output/20210928",pattern = ".csv",full.names = T)
idx <- gsub("output/20210928/","",csv_names) %>% gsub("-.*","",.) %>% as.integer()
all_idx = 1:828
table(all_idx %in% idx)
print(paste("missing",sort(as.character(all_idx[!(all_idx %in% idx)]))))
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



# Color points by dataset
# Add correlation coefficient by dataset
expr = FetchData(TASC, vars = c("SCGB1A1", "SCGB3A2","Regions"))
# Main plot
pmain <- ggplot(expr, aes_string(x = "SCGB1A1", y = "SCGB3A2",color = "Regions"))+
        geom_point()+
        scale_color_manual(values=c("#4ca64c","#E6AB02","#FF0000")) + theme_classic()
#stat_cor(aes(color = Regions), method = "spearman")

# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
        geom_density(data = expr, aes_string(x = "SCGB1A1", fill = "Regions"),
                     alpha = 0.5, size = 0.2)+
        scale_color_manual(values=c("#4ca64c","#E6AB02","#FF0000"))
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
        geom_density(data = expr, aes_string(x = "SCGB3A2", fill = "Regions"),
                     alpha = 0.5, size = 0.2)+
        coord_flip()+
        scale_color_manual(values=c("#4ca64c","#E6AB02","#FF0000"))
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

jpeg(paste0(path,"TASCs.jpeg"), units="in", width=7, height=7,res=600)
print(ggdraw(p2))
dev.off()
