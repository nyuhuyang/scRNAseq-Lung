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
