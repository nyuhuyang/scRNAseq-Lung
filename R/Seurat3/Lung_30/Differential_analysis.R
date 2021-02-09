########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
####################################
invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
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

#================== heatmap ================
object = readRDS(file = "data/Lung_30_20200710.rds")
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
object %<>% AddMetaColor(label= "annotations3", colors = Singler.colors)
Idents(object) = "annotations3"
wanted_cells = "En-L"
sub_object = subset(object, idents = wanted_cells)
Idents(sub_object) = "conditions"
read.path = "Yang/Lung_30/DE_analysis/A_Sample_types/"


deg <- read.csv(file = paste0(read.path,"Lung_30_A_312_celltypes=32_En-L_COPD_vs_distal.csv"),
                header = TRUE, row.names = 1)
deg = deg[deg$cluster == "COPD",]
Lung_markers = deg[deg$p_val_adj<0.05, "gene"]
sub_object %<>% ScaleData(features=unique(c(as.character(Lung_markers))))

DoHeatmap.1(sub_object, features = Lung_markers,do.print=T,
            angle = 0,group.bar = T, title.size = 13, no.legend = F,size=4,hjust = 0.5,
            assay = "SCT",pal_gsea = T,
            label=T, cex.row=4, legend.size = 5,width=10, height=7,unique.name = "conditions",
            title = paste("Top 160 markers of in En-L samples"))

sub_object = subset(sub_object, idents = c("distal","COPD"))

# prepare DEGs for monocle 2 ===========
# read all csv file form Yang/Lung_30/DE_analysis/A_Sample_types Epithelial section
read.path = "Yang/Lung_30/DE_analysis/A_Sample_types"

args_mt <- matrix(1:392,ncol = 7)
args_mt = args_mt[2:15,]
args_mt[args_mt <10] %<>% paste0("0",.)
args_mt %<>% paste0("Lung_30_A_",.) %>% as.vector

list_files_A <- list.files(read.path,pattern ="Lung_30_A_", full.names = T)
list_files_AE <- grep(paste(args_mt, collapse = "|"), list_files_A, value = T)

corrupted_files_A <- c()
deg_list <- list()
pb <- progress::progress_bar$new(total = length(list_files_AE))
for(i in seq_along(list_files_AE)) {
        if(file.info(list_files_AE[i])$size < 10) {
                corrupted_files_A = c(corrupted_files_A,i)
                next
        }
        df <- read.csv(list_files_AE[i],header = TRUE, row.names = 1)
        deg_list[[i]] <- df[df$p_val_adj <0.05,]
        pb$tick()
}

deg_df <- bind_rows(deg_list)
top <-  deg_df %>%
        group_by(cluster) %>%
        top_n(300, -p_val_adj)
table(top$cluster)
top = top[!duplicated(top$gene),]
table(top$cluster)
length(unique(top$gene))
write.csv(as.character(as.vector(top$gene)), file = paste0(read.path,"/top1000_epi_genes.csv"), 
          row.names = F)
# prepare DEGs for monocle 2 ===========
# read all csv file form Yang/Lung_30/DE_analysis/A_Sample_types Epithelial section
read.path = "Yang/Lung_30/DE_analysis/C_Cell_types"
args <- c(1:11,26:30,39,42,47:48)
args[args <10] %<>% paste0("0",.)

args %<>% paste0("Lung_30-",.)

list_files_C <- list.files(read.path,pattern ="Lung_30-", full.names = T)
list_files_CE <- grep(paste(args, collapse = "|"), list_files_C, value = T)

corrupted_files_C <- c()
deg_list <- list()
pb <- progress::progress_bar$new(total = length(list_files_CE))
for(i in seq_along(list_files_CE)) {
        if(file.info(list_files_CE[i])$size < 10) {
                corrupted_files_C = c(corrupted_files_C,i)
                next
        }
        df <- read.csv(list_files_CE[i],header = TRUE, row.names = 1)
        deg_list[[i]] <- df[df$p_val_adj <0.05,]
        pb$tick()
}

deg_df <- bind_rows(deg_list)
top <-  deg_df %>%
        group_by(cluster) %>%
        top_n(125, avg_logFC)
table(top$cluster)
top = top[!duplicated(top$gene),]
table(top$cluster)
length(unique(top$gene))
write.csv(as.character(as.vector(top$gene)), file = paste0(read.path,"/top1000_epi_genes.csv"), 
          row.names = F)

#============= prepare complex_cell_trajectory_ =========
read.path = "Yang/Lung_30/Monocle2/Within_regions/"
sample_string = c("distal",#128GB
                  "terminal",#32GB
                  "proximal",#32GB
                  "COPD",#32GB
                  "distal,terminal,proximal,COPD")#256GB
opts = data.frame(methods = rep(c("UseVariableGenes","ReadDE"),each = 10),
                  samples = rep(sample_string, 4),
                  root = rep(c(NA,"BC",NA,"BC"), each = 5),
                  stringsAsFactors = F)

args_string = 1:20
pb <- progress::progress_bar$new(total = length(args_string))
for(i in seq_along(args_string)){
        print(args <- opts[args_string[i],])
        Get_DE_genes <- args$methods
        sample <- stringr::str_split(args$samples, pattern = ",")[[1]]
        root = args$root
        save.path = paste0(read.path,i,"-",Get_DE_genes,'-', paste(sample, collapse = "."),"-root=",root,"/")
        cds = readRDS(paste0(save.path,basename(save.path),"_cds.rds"))
        lib_info_with_pseudo <- pData(cds)
        
        data_df <- t(monocle::reducedDimS(cds)) %>% as.data.frame() %>% 
                select_(dim_1 = 1, dim_2 = 2) %>% 
                tibble::rownames_to_column("barcodes") %>% 
                left_join(lib_info_with_pseudo %>% tibble::rownames_to_column("barcodes"), 
                          by = "barcodes")
        write.csv(data_df, file = paste0(save.path,basename(save.path),"_coordinates.csv"))
        pb$tick()
}

#================== DE by age ================
read.path = "output/20201101/"
save.path = "Yang/Lung_30/DE_analysis/D_age/"
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

csv_list <- list.files(read.path, pattern = ".csv",full.names = T)
for(i in 1:length(csv_list)){
        if(file.info(csv_list[i])$size < 10) next
        tmp = read.csv(csv_list[i], stringsAsFactors = F,row.names = 1)
        tmp = tmp[tmp$avg_logFC > 0 & tmp$p_val_adj < 0.05, ]

        write.csv(tmp, file = paste0(save.path, basename(csv_list[i])))
        svMisc::progress(i/length(csv_list)*100)
}


#======
object = readRDS(file = "data/Lung_30_20200710.rds") 
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")

Idents_list = list("older" = c("UNC-54-D", "UNC-57-D", "UNC-66-D", "UNC-70-D"),
                   "younger" = c("UNC-44-D", "UNC-48-D", "UNC-55-D", "UNC-67-D", "UNC-69-D", "UNC-71-D", "VU-27-D"))
Idents(object) = "orig.ident"
object %<>% subset(idents = unlist(Idents_list))
object@meta.data$age = "younger"
object@meta.data[object$orig.ident %in% Idents_list$older,"age"] = "older"

df_annotation <- readxl::read_excel("doc/20200903_Comparison groups of cells - Yang modified.xlsx",
                                    sheet = sub("-age","",step))
groups = df_annotation$`CELL GROUPS`
groups = groups[!is.na(groups)]
groups = gsub(" \\(.*", "", groups)
group_list <- stringr::str_split(groups, pattern = "\\+")
names(group_list) = groups
group_list[group_list == "ALL IMMUNE CELLS"][[1]] = 
        unique(unlist(group_list[35:53]))
group_list[group_list == "ALL CELLS"][[1]] = 
        unique(unlist(group_list[2:53]))
write.csv(as.data.frame.matrix(table(object$annotations3, object$age)), 
          file = paste0(save.path, "Distal_celltypes_vs_age.csv"))

# ================
save.path <- paste0("Yang/Lung_30/DE_analysis/groups/")
cell.type_list <- list("Epithelial" = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","d-S",
                                        "H","p-C","C1","C2","C3","Ion","NE","g-Muc",
                                        "g-Ser","AT1","AT2","AT2-1","AT2-p","ME"),
                       "Surface Airway Epithelial" = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","d-S",
                                                       "H","p-C","C1","C2","C3","Ion","NE"),
                       "Airway Epithelial" = c("AT1","AT2","AT2-1","AT2-p"),
                       "Structural" = c("ME","Cr","Gli","F1","F2","F3","F4",
                                        "Nr","Pr","SM1","SM2","SM3",
                                        "En-a","En-c","En-c1","En-l","En-p","En-sm","En-v"),
                       "Immune" = c("B","DC",
                                    "M-p","M0","M1","M1-2","M2","MC","Mon","Neu",
                                    "p-DC","PC","RBC",
                                    "T-cn","T-ifn","T-int","T-NK","T-p","T-reg","T-rm")
)

for(g in names(cell.type_list[c(1,3:5)])){
        csv_list = list.files(paste0(save.path,g),full.names = T)
        deg_list = pbapply::pblapply(csv_list, function(x) {
                read.csv(x,row.names = 1) %>% filter(avg_logFC > 0)
                })
        names(deg_list) = cell.type_list[[g]]
        openxlsx::write.xlsx(deg_list, 
                             file =  paste0(save.path,"DE_results_",g,".xlsx"),
                             colNames = TRUE,row.names = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))
}

