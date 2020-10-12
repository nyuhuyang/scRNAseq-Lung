########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","openxlsx",
                   "MAST","future","gplots"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

set.seed(101)
#================== DE on cluster ================
read.path = "output/20201002/"
# change the current plan to access parallelization
opts = data.frame(resolution = c(rep(0.4,16),
                                 rep(1,21),
                                 rep(2,34)),
                  cluster = c(0:15,
                              0:20,
                              0:33),
                  stringsAsFactors = F)
csv_list <- c()
for(i in 1:nrow(opts)){
    res = opts$resolution[i]
    cluster = opts$cluster[i]
    csv_list[i] <- paste0("Lung_time_6_FC1-res=",res,"_cluster=",cluster,".csv")
}
list_files <- list.files(path = read.path, pattern = "Lung_time_6_FC1-res=",full.names = F)
csv_list[!(csv_list %in% list_files)]
which(!(csv_list %in% list_files))

gde_list <- list()
resolutions = c(0.4,1,2)
for(i in seq_along(resolutions)){
    res = resolutions[i]
    clusters = opts[opts$resolution %in% res, "cluster"]
    csv_list <- paste0(read.path, "Lung_time_6_FC1-res=",res,"_cluster=",clusters,".csv")
    gde.temp <-  lapply(csv_list, read.csv) %>% bind_rows
    gde.temp = gde.temp[gde.temp$p_val_adj < 0.05,]
    gde_list[[i]] = gde.temp
}
names(gde_list) = paste0("resolutions=",c(0.4,1,2))
write.xlsx(gde_list, file = paste0(path,"Lung_time_6_DEG.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))


#================== DE on cluster ================
read.path = "output/20201006/"
# change the current plan to access parallelization
clusters =  c("BC-0","BC-m","BC-p","BC-S","BC1","BC2",
              "C1","C2","e-S","IC-S","IC1","IC2","IC2-s",
              "Ion","NEC","p-C","S","S-d","Un")
csv_list <- c()
for(i in 1:length(clusters)){
    csv_list[i] <- paste0("Lung_time_6_FC0.05_annotation=",clusters[i],".csv")
}
list_files <- list.files(path = read.path, pattern = "Lung_time_6_FC0.05_annotation=",full.names = F)
csv_list[!(csv_list %in% list_files)]
which(!(csv_list %in% list_files))

list_files <- list.files(path = read.path, pattern = "Lung_time_6_FC0.05_annotation=",
                         full.names = T)
gde_list <-  lapply(list_files, read.csv) %>% bind_rows
gde_list = gde_list[gde_list$p_val_adj < 0.05,]
gde_list = gde_list[gde_list$avg_logFC > 0 ,]

save.path = "Yang/Lung_30_time_6/Lung_time_6/DE files/"
write.xlsx(gde_list, file = paste0(save.path,"Lung_time_6_DEG_celltypes.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))


# prepare DEGs for monocle 2 ===========
# read all csv file form Yang/Lung_30/DE_analysis/A_Sample_types Epithelial section
read.path = "Yang/Lung_30_time_6/Lung_time_6/DE files/"

df_samples <- readxl::read_excel(paste0(read.path, "Lung_time_6_DEG_celltypes.xlsx"))


top <-  df_samples %>%
    group_by(cluster) %>%
    top_n(125, avg_logFC)
table(top$cluster)
top = top[!duplicated(top$gene),]
top = top[top$avg_logFC >0.7,]
length(unique(top$gene))

table(top$cluster)
range(top$avg_logFC)
write.csv(as.character(as.vector(top$gene)), file = paste0(read.path,"/top1000_epi_time6_genes.csv"), 
          row.names = F)
