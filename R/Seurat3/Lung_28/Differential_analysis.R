########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr",
                   "MAST","future","gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

read.path = "output/20200529/"
# change the current plan to access parallelization
opts = data.frame(resolution = c(rep(1,47),
                                 rep(2,64),
                                 rep(3,84),
                                 rep(4,102),
                                 rep(5,119)),
                  cluster = c(0:46,
                              0:63,
                              0:83,
                              0:101,
                              0:118),
                  stringsAsFactors = F)
set.seed(101)

csv_list <- c()
for(i in 1:nrow(opts)){
        (res = opts$resolution[i])
        (cluster = opts$cluster[i])
        csv_list[i] <- paste0("Lung_28-res=",res,"_cluster=",cluster,".csv")
}
list_files <- list.files(path = path, pattern = "Lung_28-res=")
csv_list[!(csv_list %in% list_files)]
which(!(csv_list %in% list_files))

gde_list <- list()
for(res in 1:5){
        temp_csv <- list.files(path = read.path, pattern = paste0("Lung_28-res=",res),
                                full.names = TRUE)
        l_clusters = length(temp_csv) - 1
        gde.all <- list()
        for(k in 0:l_clusters){
                gde.all[[k+1]] = read.csv(paste0(read.path,"Lung_28-res=",res,"_cluster=",k,".csv"),
                                          stringsAsFactors = F)
        }
        gde <- bind_rows(gde.all)
        gde = gde[gde$avg_logFC >0.5 & gde$p_val < 0.05,]
        gde_list[[res]] = gde
        Progress(res, 5)
}

names(gde_list) = paste0("res=",1:5)
