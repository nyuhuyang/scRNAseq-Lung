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

#==================
read.path = "output/20200703/"
# change the current plan to access parallelization
opts = data.frame(resolution = c(rep(2,75),
                                 rep(3,95),
                                 rep(4,110)),
                  cluster = c(0:74,
                              0:94,
                              0:109),
                  stringsAsFactors = F)
csv_list <- c()
for(i in 1:nrow(opts)){
        res = opts$resolution[i]
        cluster = opts$cluster[i]
        csv_list[i] <- paste0("Lung_30_FC1-res=",res,"_cluster=",cluster,".csv")
}
list_files <- list.files(path = read.path, pattern = "Lung_30_FC1-res=",full.names = F)
csv_list[!(csv_list %in% list_files)]
which(!(csv_list %in% list_files))

gde_list <- list()
resolutions = 2:4
for(r in seq_along(resolutions)){
        res = resolutions[r]
        clusters = opts[opts$resolution %in% res, "cluster"]
        csv_list <- paste0(read.path, "Lung_30_FC1-res=",res,"_cluster=",clusters,".csv")
        gde.temp <-  lapply(csv_list, read.csv) %>% bind_rows
        gde_list[[r]] = gde.temp
}
names(gde_list) = paste0("resolutions=",2:4)
write.xlsx(gde_list, file = paste0(path,"Lung_30_DEG.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

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


#==================
read.path = "output/20200624/"
opts = data.frame(resolution = c(rep(2,77),
                                 rep(4.8,119),
                                 rep(4.9,122)),
                  cluster = c(0:76,
                              0:118,
                              0:121),
                  stringsAsFactors = F)
set.seed(101)
resolutions = c(4.8, 4.9)
for(res in resolutions){
        clusters = opts[opts$resolution %in% res,"cluster"]
        csv_list <- c()
        for(i in seq_along(clusters)){
                cluster = clusters[i]
                csv_list[i] <- paste0("Lung_29_FC1-res=",res,"_cluster=",cluster,".csv")
        }
        list_files <- list.files(path = read.path, pattern = paste0("Lung_29_FC1-res=",res))
        csv_list[!(csv_list %in% list_files)] 
        print(which(!(csv_list %in% list_files))+119)

        
}




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
