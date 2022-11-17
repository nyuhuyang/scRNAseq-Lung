library(Seurat)
library(dplyr)
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

df_samples <- readxl::read_excel("doc/20220401_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()

df_samples = df_samples[grepl("P-norm|D-norm|T-norm|D-COPD", df_samples$condition),]
df_samples = df_samples[!grepl("UNC_44_P|VU_29_D|VU_35_D", df_samples$sample),]


meta.data <- readRDS("output/Lung_30_20210831_metadata_v2.rds") %>%
    select(Cell_subtype,Cell_type, UMAP_land, Family, Superfamily, Regions,orig.ident,Doublets) %>%
    dplyr::filter(Doublets == "Singlet" & Cell_subtype != "Un")

annotation <- meta.data %>% distinct(Cell_subtype, .keep_all=TRUE) %>%
    select(c("Cell_subtype","Cell_type", "UMAP_land", "Family", "Superfamily"))
    
#============= mean and sd in all cells ==========================

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

cell_perc_long$regions = plyr::mapvalues(cell_perc_long$orig.ident,
                                         from = df_samples$sample,
                                         to = df_samples$regions)
cell_perc_long$regions %<>% gsub("distal","pre-T",.)
cell_perc_long$regions %<>% factor(levels = c("proximal","pre-T","terminal","COPD"))

total_num_sd <- 
    cell_perc_long %>% 
    group_by(regions,Symbol) %>%
    summarize(mean = mean(percentage,na.rm = TRUE),
              sd = sd(percentage,na.rm = TRUE)) %>%
        ungroup() %>%
        tidyr::pivot_wider(names_from = "regions", values_from = c("mean","sd")) %>%
    tibble::column_to_rownames("Symbol") %>%
        .[df_temp$Symbol,]


write.xlsx(total_num_sd, file = paste0(path,"averages for each sample.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding")


#============= mean and sd in superFamily cells ==========================

SampleCells_Superfamily <- 
    meta.data %>% 
    dplyr::filter(Cell_subtype != "Un") %>%
    group_by(orig.ident,Superfamily) %>%
    summarize(
        SampleCells = n()
    ) %>% 
    mutate(orig.ident_Superfamily = paste0(orig.ident,"_",Superfamily)) %>%
    tibble::column_to_rownames("orig.ident_Superfamily") 

table(meta.data$orig.ident, meta.data$Superfamily) %>% as.data.frame.matrix() %>%
    tibble::rownames_to_column("orig.ident") %>% 
    pivot_longer(!orig.ident, names_to = "Superfamily",values_to = "cell.number") -> SampleCells_Superfamily



rev(c("Cell_subtype","Cell_type", "UMAP_land", "Family", "Superfamily")) %>%
    lapply(function(type) {
        table(meta.data[,type],meta.data$orig.ident) %>% as.data.frame.matrix() %>%
            tibble::rownames_to_column("Symbol")
    }) %>% bind_rows() %>% 
    distinct(Symbol, .keep_all=TRUE) %>%
    #arrange_(factor(Symbol,levels = df_temp$Symbol)) %>%
    tibble::column_to_rownames("Symbol") %>%
    .[df_temp$Symbol,]-> cell_number


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

cell_perc_long$regions = plyr::mapvalues(cell_perc_long$orig.ident,
                                         from = df_samples$sample,
                                         to = df_samples$regions)
cell_perc_long$regions %<>% gsub("distal","pre-T",.)
cell_perc_long$regions %<>% factor(levels = c("proximal","pre-T","terminal","COPD"))

total_num_sd <- 
    cell_perc_long %>% 
    group_by(regions,Symbol) %>%
    summarize(mean = mean(percentage,na.rm = TRUE),
              sd = sd(percentage,na.rm = TRUE)) %>%
    ungroup() %>%
    tidyr::pivot_wider(names_from = "regions", values_from = c("mean","sd")) %>%
    tibble::column_to_rownames("Symbol") %>%
    .[df_temp$Symbol,]
table(rownames(total_num_sd) %in% df_temp$Symbol)
rownames(total_num_sd)[!(rownames(total_num_sd) %in% df_temp$Symbol)]

write.xlsx(total_num_sd, file = paste0(path,"averages for each sample_by_Superfamily.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding")

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
cell_perc_long$regions = plyr::mapvalues(cell_perc_long$orig.ident,
                                         from = df_samples$sample,
                                         to = df_samples$regions)
cell_perc_long$regions %<>% gsub("distal","pre-T",.)
cell_perc_long$regions %<>% factor(levels = c("proximal","pre-T","terminal","COPD"))

total_num_sd <- 
    cell_perc_long %>% 
    group_by(regions,Symbol) %>%
    summarize(mean = mean(percentage,na.rm = TRUE),
              sd = sd(percentage,na.rm = TRUE)) %>%
    ungroup() %>%
    tidyr::pivot_wider(names_from = "regions", values_from = c("mean","sd")) %>%
    tibble::column_to_rownames("Symbol") %>%
    .[df_temp$Symbol,]
table(rownames(total_num_sd) %in% df_temp$Symbol)
rownames(total_num_sd)[!(rownames(total_num_sd) %in% df_temp$Symbol)]

write.xlsx(total_num_sd, file = paste0(path,"averages for each sample_by_Family.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding")


#-======================
#after obtain DEGs from human lung HLCA v2
library(arrow)
annotations.path <- "data/human_lung_v2/data/annotations.parquet"
annotations <- read_parquet(annotations.path)
Secretory <- unique(annotations[annotations$ann_level_3 %in% c("Secretory"),"ann_finest_level"])

#-----ven diagram-------
DEGsv2 <- read_excel("output/20221013/2022-10-13-ann_finest_level.xlsx") %>%
    filter(p_val_adj < 0.05)
colnames(DEGsv2) =c("cluster","gene","scores","avg_log2FC","p_val","p_val_adj","pct.1","pct.2")
DEGsv2 =DEGsv2[,c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene")]
openxlsx::write.xlsx(DEGsv2, file =  "output/20221013/2022-10-13-ann_finest_level.filter.xlsx",
                     colNames = TRUE,rowNames = F,borders = "surrounding")

DEGsv2 <- read_excel("output/20221013/2022-10-13-ann_finest_level.xlsx") %>% 
    filter(p_val_adj < 0.05 & pts > 0.1 & avg_log2FC >1 & group %in% Secretory)  %>%
    group_by(group) %>%
    top_n(100,scores) %>%
    ungroup() %>%
    as.data.frame()
Secretory_pos_genes <- split(DEGsv2,DEGsv2$group) %>% 
    lapply(function(df) df[(df$avg_log2FC > 1),"genes"]) %>% 
    lapply(unique)


DEGs <- read_excel("Yang/Lung_30/hg38/DE_analysis/Lung_30_DEG_Cell.category.xlsx") %>% as.data.frame()
DEGs[DEGs$p_val == 0,"p_val"] = min(DEGs[DEGs$p_val > 0,"p_val"], .Machine$double.xmin)
pos_genes <-  DEGs %>%
    filter(p_val_adj < 0.05 & cluster == "TASC") %>% 
    head(100) %>%
    .[,"gene"] %>% unique
Secretory_pos_genes %<>% c(list("TASC"=pos_genes))

for(cell in Secretory){
    euler_df <- eulerr::euler(Secretory_pos_genes[c("TASC",cell)],shape = "circle")
    
    g <- plot(euler_df, quantities = TRUE, lty = 1:6,
              legend = TRUE)
    jpeg(paste0(path,"/","Venn_",cell,"_TASC_top100.jpeg"), units="in", width=10, height=7,res=600)
    print(g)
    dev.off()
    Progress(which(Secretory %in% cell), length(Secretory))
}

library(ggvenn)
ggvenn(
    Secretory_pos_genes, 
    #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 8
)
x <- list(
    A = sample(genes,300), 
    B = sample(genes,525), 
    C = sample(genes,440),
    D = sample(genes,350),
    E = sample(genes,250),
    EF = sample(genes,250)
)
venn.diagram(Secretory_pos_genes[1,], filename = paste0(path,"venn_8.png"))

library(VennDiagram)
openxlsx::write.xlsx(DEGsv2, file =  "output/20221013/2022-10-13-ann_finest_level.Secretory.xlsx",
                     colNames = TRUE,rowNames = F,borders = "surrounding")

DEGs_AT2 <- read_excel("output/20221014/2022-10-14-Transitional Club-AT2 vs AT2 in   .xlsx") %>% as.data.frame() %>%
    filter(p_val_adj < 0.05) %>% 
    head(100) %>%
    .[,"gene"] %>% unique
DEGs_Secretory <- read_excel("output/20221014/2022-10-14-Transitional Club-AT2 vs Club (nasal)_Club (non-nasal)_Goblet (bronchial)_Goblet (nasal)_Goblet (subsegmental) in   .xlsx") %>% as.data.frame() %>%
    filter(p_val_adj < 0.05) %>% 
    head(100) %>%
    .[,"gene"] %>% unique

Secretory_AT2_pos_genes <- list("TASC"=pos_genes,
                                "Transitional Club-AT2 vs AT2" = DEGs_AT2,
                                "Transitional Club-AT2 vs other Secretory" = DEGs_Secretory)
euler_df <- eulerr::euler(Secretory_AT2_pos_genes,shape = "circle")

g <- plot(euler_df, quantities = TRUE, lty = 1:6,
          legend = TRUE)
file.name = "Venn_Secretory_AT2_TASC_top100.jpeg"
jpeg(paste0(path,"/",file.name), units="in", width=10, height=7,res=600)
print(g)
dev.off()



#-----fgsea dotplot-------
# all
DEGsv2 <- read_excel("output/20221013/2022-10-13-ann_finest_level.xlsx") %>% 
    filter(p_val_adj < 0.05 & pts > 0.1 & avg_log2FC >2)  %>%
    group_by(group) %>%
    top_n(250) %>%
    ungroup()
    
pos_genes <- split(DEGsv2,DEGsv2$group) %>% 
    lapply(function(df) df$genes) %>% 
    lapply(unique)
str(pos_genes)
res <- read_excel("Yang/Lung_30/hg38/DE_analysis/Lung_30_DEG_Cell.category.xlsx") %>% 
    filter(pct.1 > 0.1 & avg_log2FC >0) %>% 
    as.data.frame() 
plot <- FgseaDotPlot(stats=res, pathways=pos_genes,Rowv = TRUE,Colv = TRUE,
             cols = c(rev(brewer.pal(10,"RdBu"))), return.plot = TRUE)
plot <- plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
             axis.text.x =    element_text(angle = 90,
                                           hjust = 0,
                                           vjust =  0.5)
             )
jpeg(paste0(path, "/Dotplot_log2FC_2_fdr0.25_p0.05.jpeg"),units="in", width=9, height=9,res=600)
print(plot)
dev.off()


#---
annotations.path <- "data/human_lung_v2/data/annotations.parquet"
annotations <- read_parquet(annotations.path)
Epi <- unique(annotations[annotations$ann_level_1 %in% "Epithelial","ann_finest_level"])
Epi_pos_genes <- pos_genes[Epi]

meta.data = readRDS("output/Lung_30_20210831_metadata_v2.rds")
Lung_Epi <- unique(meta.data[meta.data$Superfamily %in% "Epi","Cell_subtype"]) %>% as.character()

plot1 <- FgseaDotPlot(stats= res[res$cluster %in% Lung_Epi,], pathways=Epi_pos_genes,Rowv = TRUE,Colv = TRUE,cols = rev(brewer.pal(10,"RdBu")), return.plot = TRUE,font.main = 20)

plot1 <- plot1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                     axis.text.x =    element_text(angle = 90,
                                                   hjust = 0,
                                                   vjust =  0.5)
)
jpeg(paste0(path, "/Dotplot_Epi_log2FC_2_fdr0.25_p0.05.jpeg"),units="in", width=5, height=5 ,res=600)
print(plot1)
dev.off()
