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

read.path = "output/20200619/"
# change the current plan to access parallelization
opts = data.frame(cluster = c(0:76),
                  stringsAsFactors = F)
set.seed(101)
res = 2
csv_list <- c()
for(i in 1:nrow(opts)){
        (cluster = opts$cluster[i])
        csv_list[i] <- paste0("Lung_29-res=",res,"_cluster=",cluster,".csv")
}
list_files <- list.files(path = read.path, pattern = "Lung_29-res=")
csv_list[!(csv_list %in% list_files)]
which(!(csv_list %in% list_files))

gde_list <- list()
temp_csv <- list.files(path = read.path, pattern = paste0("Lung_29-res=",res),
                       full.names = TRUE)
l_clusters = length(temp_csv) - 1
gde.all <- list()
for(k in 0:l_clusters){
        gde.all[[k+1]] = read.csv(paste0(read.path,"Lung_29-res=",res,"_cluster=",k,".csv"),
                                  stringsAsFactors = F)
        Progress(k, l_clusters+1)
}
gde <- bind_rows(gde.all)
gde = gde[gde$avg_logFC >0.5 & gde$p_val < 0.05,]
gde_list[[1]] = gde
        

names(gde_list) = paste0("res=",2)
write.xlsx(gde_list, file = paste0(path,"DEG_markers_Nointegration.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
