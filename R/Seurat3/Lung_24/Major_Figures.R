########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","tidyr","magrittr","ggpubr","MAST"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#######################
# box plot
#######################

# for major cell ===========================
# load data
cell_dist <- readxl::read_excel("doc/Cell distribution - for plots.xlsx")
cell_dist = cell_dist[!is.na(cell_dist$Proximal),]
colnames(cell_dist) = gsub(" ",".",cell_dist[1,])
cell_dist = cell_dist[-1,]
cell_dist = cell_dist[,-grep("Families|Code",colnames(cell_dist))]

regions = c("Proximal","Distal","Terminal")
major_cells <- c("Epithelial","Stromal","Endothelial","Immune","Neuronal")
major_cell_dist = as.data.frame(cell_dist[1:5,])
df_major_cell_dist <- gather(major_cell_dist, key = "samples", value = "percentage",-Cell.categories)
head(df_major_cell_dist)
df_major_cell_dist$percentage %<>% as.numeric() 
df_major_cell_dist$regions = gsub("-R","",df_major_cell_dist$samples)
df_major_cell_dist$regions %<>% gsub(".*-","",.)
df_major_cell_dist$regions %<>% plyr::mapvalues(from = c("P","D","T"),
                                       to = regions)
df_major_cell_dist$regions %<>% as.factor()
df_major_cell_dist$regions %<>% factor(levels = rev(regions))
table(df_major_cell_dist$regions)
#df_major_cell_dist = df_major_cell_dist[,-grep("samples",colnames(df_major_cell_dist))]
#rownames(df_major_cell_dist) = df_major_cell_dist$Cell.categories
df_major_cell_dist1 <- spread(df_major_cell_dist, Cell.categories,percentage, fill = 0)
head(df_major_cell_dist1)

jpeg(paste0(path,"Cell distribution.jpeg"), units="in", width=10, height=7,res=600)
ggboxplot(df_major_cell_dist1, x = "regions",
          y = rev(major_cells),
          merge = "flip",
          ylab = "Cell percenrage", 
          xlab = "",
          palette = c('#3D9970','#FF4136','#FF851B'),
          color = "regions",
          add = "jitter",
          #add.params = list(size = 0.1, jitter = 0.2),
          rotate = TRUE)
dev.off()

my_comparisons <- list(c("Proximal", "Distal"), c("Distal", "Terminal"),c("Proximal", "Terminal"))

for(c in major_cells){
        g <- ggboxplot(df_major_cell_dist1, x = "regions",
                       y = c,
                       merge = "flip",
                       ylab = "Cell percenrage", 
                       xlab = "",
                       palette = c('#3D9970','#FF4136','#FF851B'),
                       color = "regions",
                       add = "jitter",
                       #add.params = list(size = 0.1, jitter = 0.2),
                       rotate = TRUE)+
                stat_compare_means(comparisons = my_comparisons)
        jpeg(paste0(path,"Cell distribution-",c,".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
}

# for minor cell =======================
# load data
cell_dist <- readxl::read_excel("doc/Cell distribution - for plots.xlsx")
cell_dist = cell_dist[!is.na(cell_dist$Proximal),]
colnames(cell_dist) = gsub(" ",".",cell_dist[1,])
cell_dist = cell_dist[-1,]
cell_dist = cell_dist[,-grep("Cell.categories",colnames(cell_dist))]
#https://stackoverflow.com/questions/10554741/fill-in-data-frame-with-values-from-rows-above/32536507
cell_dist$Families = zoo::na.locf(cell_dist$Families)
cell_dist = as.data.frame(cell_dist[-c(1:5),])

regions = c("Proximal","Distal","Terminal")
(minor_cells <- unique(cell_dist$Families))
for(cell in minor_cells) {
        print(cell)
        cell_path <- paste0(path,cell,"/")
        if(!dir.exists(cell_path))dir.create(cell_path, recursive = T)
        
        minor_cell_dist = cell_dist[cell_dist$Families %in% cell,-1]
        df_minor_cell_dist <- gather(minor_cell_dist, key = "samples", value = "percentage",-Code)
        head(df_minor_cell_dist,2)
        df_minor_cell_dist$percentage %<>% as.numeric() 
        df_minor_cell_dist$regions = gsub("-R","",df_minor_cell_dist$samples)
        df_minor_cell_dist$regions %<>% gsub(".*-","",.)
        df_minor_cell_dist$regions %<>% plyr::mapvalues(from = c("P","D","T"),
                                                        to = regions)
        df_minor_cell_dist$regions %<>% as.factor()
        df_minor_cell_dist$regions %<>% factor(levels = rev(regions))
        table(df_minor_cell_dist$regions)
        df_minor_cell_dist1 <- spread(df_minor_cell_dist, Code,percentage, fill = 0)
        ggboxplot(df_minor_cell_dist1, x = "regions",
                  y = rev(minor_cell_dist$Code),
                  bxp.errorbar = TRUE,
                  merge = "flip",
                  width = 0.7,
                  size = NULL,
                  ylab = "Cell percenrage", 
                  xlab = "",
                  palette = c('#3D9970','#FF4136','#FF851B'),
                  color = "regions",
                  add = "jitter",
                  #add.params = list(size = 0, jitter = 0),
                  rotate = TRUE)
        jpeg(paste0(cell_path,"Cell distribution-",cell,".jpeg"), units="in", width=10, height=20,res=600)
        
        print(g)
        dev.off()
        
        my_comparisons <- list(c("Proximal", "Distal"), c("Distal", "Terminal"),c("Proximal", "Terminal"))
        colnames(df_minor_cell_dist1) %<>% gsub("-|/","_",.)
        (minor_cell_types <- colnames(df_minor_cell_dist1)[-c(1:2)])

        for(c in minor_cell_types){
                df_minor_cell_dist2 <- df_minor_cell_dist1[,c("samples", "regions",c)]
                g1 <- ggboxplot(df_minor_cell_dist1, x = "regions",
                               y = c,
                               merge = "flip",
                               ylab = "Cell percenrage", 
                               xlab = "",
                               palette = c('#3D9970','#FF4136','#FF851B'),
                               color = "regions",
                               add = "jitter",
                               #add.params = list(size = 0.1, jitter = 0.2),
                               rotate = TRUE)+
                        stat_compare_means(comparisons = my_comparisons)
                jpeg(paste0(cell_path,"Cell distribution-",c,".jpeg"), units="in", width=10, height=7,res=600)
                print(g1)
                dev.off()
        }
        Progress(which(minor_cells %in% cell),length(minor_cells))
}

#######################
# GSEA
#######################
csv_files <- list.files(path,pattern = "Lung_24-FC0_markers_")

        paste0(path,"Lung_24-FC0_markers_",args,"_",cell.type,".csv"))
