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

#To address minor comment 3 of Reviewer 1, we need error bars (standard deviation) for % of individual cell subtypes (% of all cells per sample) in normal pre-T and COPD pre-T samples.
#Attached is the table that contains the original data which you generated. I highlighted the cell subtypes and columns containing the mean (average) values (% all cell types) for pre-T and COPD pre-T groups. We need standard deviation data for these.
df_temp <- readxl::read_excel("Yang/TASC paper/STab-03 - TASC - S Table III. Cell distribution.xlsx",
                              sheet = "Sheet1") %>%
    filter(Symbol != "Un")

df_samples <- readxl::read_excel("doc/202108014_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples$smoking %<>% gsub(".*PY|Cigars", "S",.)
df_samples = df_samples[grepl("D-norm", df_samples$condition),]
df_samples = df_samples[!grepl("VU_29_D|VU_35_D", df_samples$sample),]

# load data
meta.data <- readRDS("output/Lung_30_20210831_metadata_v2.rds") %>%   
    subset(select=which(!duplicated(names(.)))) %>%
    dplyr::filter(Doublets == "Singlet" & Cell_subtype != "Un"  & smoking == "pre-T") %>%
    select(Cell_subtype,Cell_type, UMAP_land, Family, Superfamily, smoking,orig.ident,Doublets)

meta.data$orig.ident %<>% droplevels()
meta.data$smoking = plyr::mapvalues(meta.data$orig.ident,
                                      from = df_samples$sample,
                                      to = df_samples$smoking)

annotation <- meta.data %>% distinct(Cell_subtype, .keep_all=TRUE) %>%
    select(c("Cell_subtype","Cell_type", "UMAP_land", "Family", "Superfamily"))

#============= mean in all cells ==========================

SampleCells <- meta.data %>% 
    dplyr::filter(Cell_subtype != "Un") %>%
    group_by(orig.ident) %>%
    summarize(
        SampleCells = n()
    ) %>% tibble::column_to_rownames("orig.ident") %>% t

rev(c("Cell_subtype","Cell_type", "UMAP_land", "Family", "Superfamily")) %>%
    lapply(function(type) {
        table(meta.data[,type],meta.data$orig.ident) %>% as.data.frame.matrix() %>%
            tibble::rownames_to_column("Symbol")
    }) %>% bind_rows() %>% 
    distinct(Symbol, .keep_all=TRUE) %>%
    #arrange_(factor(Symbol,levels = df_temp$Symbol)) %>%
    tibble::column_to_rownames("Symbol") %>%
    .[df_temp$Symbol,]-> cell_number

cell_perc_total <- sweep(cell_number, 2, SampleCells,"/") *100
cell_perc_total %>% tibble::rownames_to_column("Symbol") %>% 
    pivot_longer(!Symbol, names_to = "orig.ident", values_to = "percentage") ->
    cell_perc_long

cell_perc_long$smoking = plyr::mapvalues(cell_perc_long$orig.ident,
                                         from = df_samples$sample,
                                         to = df_samples$smoking)

total_num_sd <- 
    cell_perc_long %>% 
    group_by(smoking,Symbol) %>%
    summarize(mean = mean(percentage,na.rm = TRUE),
              sd = sd(percentage,na.rm = TRUE)) %>%
    ungroup() %>%
    tidyr::pivot_wider(names_from = "smoking", values_from = c("mean","sd")) %>%
    tibble::column_to_rownames("Symbol") %>%
    .[df_temp$Symbol,]
write.xlsx(total_num_sd, file = paste0(path,"averages for each sample.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding")
#===========
S_NS_wilcox.test <- function(df = cell_perc_long, group1 = "S",
                               group2 = "NS"){
    
    df = df[df$smoking %in% unique(c(group1,group2)),]
    df1 = df %>% dplyr::filter(smoking %in% group1) %>% pivot_wider(id_cols = !c("smoking"), 
                                                                    names_from = "orig.ident", 
                                                                    values_from = "percentage") %>% 
        column_to_rownames(var="Symbol")
    df2 = df %>% dplyr::filter(smoking %in% group2) %>% pivot_wider(id_cols = !c("smoking"), 
                                                                    names_from = "orig.ident", 
                                                                    values_from = "percentage") %>% 
        column_to_rownames(var="Symbol")
    sdf1 <- split(df1, rownames(df1)) %>% 
        rapply(f=function(x) ifelse(is.nan(x),0,x), how="replace" )
    sdf2 <- split(df2, rownames(df2)) %>% 
        rapply(f=function(x) ifelse(is.nan(x),0,x), how="replace" )
    
    p_value = mapply(function(x,y) {
        scores = wilcox.test(x= unlist(x), y= unlist(y), alternative = "two.sided",exact=FALSE,
                             na.action = na.omit)
        return(scores$p.value)
    }, sdf1, sdf2)
    
    return(p_value)
}

wilcox_res_all = S_NS_wilcox.test(cell_perc_long, group1 = "S",group2 = "NS")


#============= mean and sd in superFamily cells ==========================
table(meta.data$orig.ident, meta.data$Superfamily) %>% as.data.frame.matrix() %>%
    tibble::rownames_to_column("orig.ident") %>% 
    pivot_longer(!orig.ident, names_to = "Superfamily",values_to = "cell.number") -> SampleCells_Superfamily


c("Epi","Im","Str") %>%
    sapply(function(type){
        annotation %>% filter(Superfamily == type) %>% 
            pivot_longer(everything(), values_to = "Symbol") %>%
            distinct(Symbol, .keep_all=TRUE) %>%
            .[,"Symbol"] %>% pull -> cell_name
        if( any(c("BC-S","SM-Pr") %in% cell_name) )
            cell_name = cell_name[-which(cell_name %in% c("BC-S","SM-Pr"))] 
        cell_name
        # Remove multiple list elements
    }) -> Superfamily_list

table(unlist(Superfamily_list) %in% rownames(cell_number))
unlist(Superfamily_list)[!unlist(Superfamily_list) %in% rownames(cell_number)]

c("Epi","Im","Str") %>%
    lapply(function(type){
        SampleCells_Superfamily %>% filter(Superfamily == type) %>%
            select(-Superfamily)  %>% 
            #remove_rownames %>%
            tibble::column_to_rownames("orig.ident") %>% t %>%
            .[,colnames(cell_number)] ->  sub_SampleCells_Superfamily
        cell_number[Superfamily_list[[type]],] %>%
            sweep(2, sub_SampleCells_Superfamily,"/") *100
    }) %>% bind_rows() %>%
    .[df_temp$Symbol,] %>% 
    tibble::rownames_to_column("Symbol") %>% 
    pivot_longer(!Symbol, names_to = "orig.ident", values_to = "percentage") ->
    cell_perc_long

cell_perc_long$smoking = plyr::mapvalues(cell_perc_long$orig.ident,
                                         from = df_samples$sample,
                                         to = df_samples$smoking)

total_num_sd <- 
    cell_perc_long %>% 
    group_by(smoking,Symbol) %>%
    summarize(mean = mean(percentage,na.rm = TRUE),
              sd = sd(percentage,na.rm = TRUE)) %>%
    ungroup() %>%
    tidyr::pivot_wider(names_from = "smoking", values_from = c("mean","sd")) %>%
    tibble::column_to_rownames("Symbol") %>%
    .[df_temp$Symbol,]
table(rownames(total_num_sd) %in% df_temp$Symbol)
rownames(total_num_sd)[!(rownames(total_num_sd) %in% df_temp$Symbol)]

write.xlsx(total_num_sd, file = paste0(path,"averages for each sample_by_Superfamily.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding")

wilcox_res_SuperFamily = S_NS_wilcox.test(cell_perc_long, group1 = "S",group2 = "NS")

#============= mean and sd in Family cells ==========================
table(meta.data$orig.ident, meta.data$Family) %>% as.data.frame.matrix() %>%
    tibble::rownames_to_column("orig.ident") %>% 
    pivot_longer(!orig.ident, names_to = "Family",values_to = "cell.number") -> SampleCells_Family

rev(c("Cell_subtype","Cell_type", "UMAP_land", "Family", "Superfamily")) %>%
    lapply(function(type) {
        table(meta.data[,type],meta.data$orig.ident) %>% as.data.frame.matrix() %>%
            tibble::rownames_to_column("Symbol")
    })  %>% bind_rows() %>% 
    distinct(Symbol, .keep_all=TRUE) %>%
    #arrange_(factor(Symbol,levels = df_temp$Symbol)) %>%
    tibble::column_to_rownames("Symbol") %>%
    .[df_temp$Symbol,]-> cell_number


c("ASE","En","Im-l","Im-m","SMG","Strm","AT" ) %>%
    sapply(function(type){
        annotation %>% filter(Family == type) %>% 
            pivot_longer(everything(), values_to = "Symbol") %>%
            distinct(Symbol, .keep_all=TRUE) %>%
            .[,"Symbol"] %>% pull -> cell_name
        if( any(c("BC-S","SM-Pr","Epi","Str","Im") %in% cell_name) )
            cell_name = cell_name[-which(cell_name %in% c("BC-S","SM-Pr","Epi","Str","Im"))] 
        cell_name
        # Remove multiple list elements
    }) -> Family_list

table(unlist(Family_list) %in% rownames(cell_number))
unlist(Family_list)[!unlist(Family_list) %in% rownames(cell_number)]

c("ASE","En","Im-l","Im-m","SMG","Strm","AT" ) %>%
    lapply(function(type){
        SampleCells_Family %>% filter(Family == type) %>%
            select(-Family)  %>% 
            tibble::column_to_rownames("orig.ident")  %>% t   ->  sub_SampleCells_Family
        cell_number[Family_list[[type]],] %>%
            sweep(2, sub_SampleCells_Family[1,],"/") *100
    }) %>% bind_rows() %>%
    .[df_temp$Symbol,] %>% 
    tibble::rownames_to_column("Symbol") %>% 
    pivot_longer(!Symbol, names_to = "orig.ident", values_to = "percentage") ->
    cell_perc_long
cell_perc_long <- cell_perc_long[!is.na(cell_perc_long$percentage),]
cell_perc_long$smoking = plyr::mapvalues(cell_perc_long$orig.ident,
                                         from = df_samples$sample,
                                         to = df_samples$smoking)
total_num_sd <- 
    cell_perc_long %>% 
    group_by(smoking,Symbol) %>%
    summarize(mean = mean(percentage,na.rm = TRUE),
              sd = sd(percentage,na.rm = TRUE)) %>%
    ungroup() %>%
    tidyr::pivot_wider(names_from = "smoking", values_from = c("mean","sd")) %>%
    tibble::column_to_rownames("Symbol") %>%
    .[df_temp$Symbol,]
table(rownames(total_num_sd) %in% df_temp$Symbol)
rownames(total_num_sd)[!(rownames(total_num_sd) %in% df_temp$Symbol)]

write.xlsx(total_num_sd, file = paste0(path,"averages for each sample_by_Family.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding")

wilcox_res_Family = S_NS_wilcox.test(cell_perc_long, group1 = "S",group2 = "NS")

names(wilcox_res_all)[!(names(wilcox_res_all) %in% names(wilcox_res_Family))]
wilcox_res_Family %<>% c(c("Epi"=NaN,"Im" = NaN,  "Str"=NaN))
wilcox_res_Family <- wilcox_res_Family[names(wilcox_res_SuperFamily)]

table(names(wilcox_res_all) == names(wilcox_res_SuperFamily))
table(names(wilcox_res_Family) == names(wilcox_res_SuperFamily))

wilcox_res <- bind_cols(list(wilcox_res_all,wilcox_res_SuperFamily,wilcox_res_Family)) %>% as.data.frame()
colnames(wilcox_res) = c("all","SuperFamily","Family")
rownames(wilcox_res) = names(wilcox_res_all)
wilcox_res = wilcox_res[df_temp$Symbol,]
write.xlsx(wilcox_res, file = paste0(path,"wilcox_res.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding")
