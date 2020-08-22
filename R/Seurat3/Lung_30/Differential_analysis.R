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


#================== DE on cell types ================
read.path = "output/20200814/"
args=1:62
args[args < 10] = paste0("0", args[args < 10])
cell_types = sort(c("AT1","AT2","AT2-1","AT2-p","BC","BC-p","BC-S","IC1","IC2","IC-S","H","p-C",
                    "C1","C2","C3","S","S-d","Ion","NEC","SMG-Muc","SMG-Ser","MEC","Cr",
                    "F1","F2","F3","F4","Gli","SM1","SM2","SM3","Pr","En-A","En-C","En-C1",
                    "En-V","En-p","En-SM","En-L","Nr","Neu","MC","Mon","M0","M1","M2",
                    "M1-2","M-p","DC","P-DC","B","PC","T-cn","T-reg","T-rm","T-NK","T7",
                    "T-ifn","T-int","T-p","T-un","RBC"))
csv_list <- paste0("Lung_30-",args,"_FC0.1_",cell_types,".csv")

list_files <- list.files(path = read.path, pattern ="Lung_30-")
csv_list[!(csv_list %in% list_files)]

gde.all <- lapply(paste0(read.path,list_files), function(x) {
        tmp = read.csv(x, stringsAsFactors = F)
        tmp = tmp[tmp$p_val <0.05, ]
        tmp[,-1]
        })
names(gde.all) = cell_types
gde <- bind_rows(gde.all)
write.xlsx(gde, file = paste0(path,"DEG_markers_by_cell_types.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#================== DE on group -A ================
read.path = "Yang/Lung_30/DE_analysis/"
list_files <- list.files(path = paste0(read.path,"A_Sample_types"),
                         pattern ="Lung_30_A_")
int <- gsub("Lung_30_A_","",list_files) %>% gsub("_.*","",.) %>% as.integer()
l <- 1:378
l[!(l %in% int)]

list_files_A <- list.files(path = paste0(read.path,"A_Sample_types"),
                         pattern ="Lung_30_A_",full.names = T)

corrupted_files_A <- c()
for(i in seq_along(list_files_A)) {
        if(file.info(list_files_A[i])$size < 10) {
                corrupted_files_A = c(corrupted_files_A,i)
                next
        }
        df <- read.csv(list_files_A[i],header = TRUE, row.names = 1)
        df <- df[df$p_val_adj <0.05,]
        write.csv(df, file = list_files_A[i])
        Progress(i, length(list_files_A))
}
#================== DE on group -B ================
read.path = "Yang/Lung_30/DE_analysis/"
list_files <- list.files(path = paste0(read.path,"B_Cell_groups"),
                         pattern ="Lung_30_B_")
int <- gsub("Lung_30_B_","",list_files) %>% gsub("_.*","",.) %>% as.integer()
table(1:94 %in% int)

list_files_B <- list.files(path = paste0(read.path,"B_Cell_groups"),
                           pattern ="Lung_30_B_",full.names = T)

corrupted_files_B <- c()
for(i in seq_along(list_files_B)) {
        if(file.info(list_files_B[i])$size < 10) {
                corrupted_files_B = c(corrupted_files_B,i)
                next
        }
        df <- read.csv(list_files_B[i],header = TRUE, row.names = 1)
        df <- df[df$p_val_adj <0.05,]
        write.csv(df, file = list_files_B[i])
        Progress(i, length(list_files_B))
}

#================== DE on group - C_Cell_types ================
read.path = "Yang/Lung_30/DE_analysis/"
list_files <- list.files(path = paste0(read.path,"C_Cell_types"),
                         pattern ="Lung_30-")
int <- gsub("Lung_30-","",list_files) %>% gsub("_.*","",.) %>% as.integer()
table(1:62 %in% int)

list_files_C <- list.files(path = paste0(read.path,"C_Cell_types"),
                           pattern ="Lung_30-",full.names = T)

corrupted_files_C <- c()
for(i in seq_along(list_files_C)) {
        if(file.info(list_files_C[i])$size < 10) {
                corrupted_files_C = c(corrupted_files_C,i)
                next
        }
        df <- read.csv(list_files_C[i],header = TRUE, row.names = 1)
        df <- df[df$p_val_adj <0.05,]
        write.csv(df, file = list_files_C[i])
        Progress(i, length(list_files_C))
}
save(list_files_A, list_files_B, list_files_C,corrupted_files_A,file= paste0("output/20200819/DE_list_files",".Rda"))
