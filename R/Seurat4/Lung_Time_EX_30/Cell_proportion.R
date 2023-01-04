library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(stringr)
library(openxlsx)
library(tidyverse)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)


# load data
meta.data = readRDS("output/20211020/meta.data_Cell_subtype_time6.rds")
table(meta.data$Doublets)
meta.data %<>% subset(Doublets %in% "Singlet")
meta.data$day %<>% factor(levels = paste0("D",c(0,3,7,14,21,28)))

#============ B - Cell numbers in each sample ============================
Cell_types <- c("Cell_subtype","Cell_type","Major_Cell_type")
cell_number = lapply(Cell_types, function(Cell_type){
    as.data.frame.matrix(table(meta.data[,Cell_type], meta.data$day))
}) %>% bind_rows()

Symbol = gsub("\\..*","",rownames(cell_number))
cell_number = cell_number[!duplicated(Symbol),]
rownames(cell_number) %<>% gsub("\\..*","",.)

#====
df_number <- readxl::read_excel("doc/Table II. Cell ontology and distribution - template.xlsx",
                                    sheet = "B - Cell numbers in each sample")
colnames(df_number) = df_number[2,]
df_number = df_number[-2,]
df_number = df_number[!is.na(df_number$Symbol),]
Symbol = df_number$Symbol = gsub(" .*","",df_number$Symbol)
setdiff(rownames(cell_number), Symbol) #"TNK"   "BC-S"  "SM-Pr"
setdiff(Symbol,rownames(cell_number)) #

day =  colnames(cell_number)
rownames(df_number) = Symbol
table(colnames(df_number)[4:ncol(df_number)] == colnames(cell_number))

write.xlsx(cell_number[Symbol,], file = paste0("output/20211020/","B_Cell_numbers.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

day_number = as.data.frame(table(meta.data$day))
colnames(day_number) = c("day","total")
#============ A - Cell proportions per group ============================
df_number <- readxl::read_excel("doc/Table II. Cell ontology and distribution - template.xlsx",
                                sheet = "A - Cell proportions per group")
colnames(df_number) = df_number[2,]
df_number = df_number[-2,]
df_number = df_number[!is.na(df_number$Symbol),!is.na(colnames(df_number))]
Symbol = df_number$Symbol = gsub(" .*","",df_number$Symbol)

cell_number = readxl::read_excel(path = paste0("output/20211020/","B_Cell_numbers.xlsx"),)
colnames(cell_number)[1] = "cell_types"
cell_number_long = pivot_longer(cell_number, !cell_types,names_to = "day", values_to = "num")
cell_number_long = cell_number_long[,c("cell_types","day","num")]
# mean_res % among all cells (per sample)

cell_number_long$total = plyr::mapvalues(cell_number_long$day,
                                    from = day_number$day,
                                    to = day_number$total)
cell_number_long$total %<>% as.numeric()
cell_number_long$mean = cell_number_long$num / cell_number_long$total *100
mean_res_all = cell_number_long %>% pivot_wider(!c("num","total"), names_from = "day", values_from = "mean")
mean_res_all %<>% column_to_rownames(var="cell_types")


# wilcox_res % among all cells (per sample)	
cell_number %<>% tibble::column_to_rownames("cell_types")

sdf1 <- split(cell_number[,c("D0","D3","D7")], Symbol)
sdf2 <- split(cell_number[,c("D14","D21","D28")], Symbol)

p_value = mapply(function(x,y) {
    scores = wilcox.test(x= unlist(x), y= unlist(y), alternative = "two.sided",exact=FALSE,
                         na.action = na.omit)
    return(scores$p.value)
}, sdf1, sdf2)

mean_res_all$p_value = p_value[Symbol]
 write.xlsx(mean_res_all, file = paste0(path,"A_Cell_proportions.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
