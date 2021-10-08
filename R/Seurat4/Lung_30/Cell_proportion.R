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
object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")
table(rownames(object@meta.data) == rownames(meta.data))
table(object$barcode ==meta.data$barcode)
object@meta.data = meta.data

DefaultAssay(object) = "SCT"
object %<>% subset(subset = Doublets %in% "Singlet")

df_samples <- readxl::read_excel("doc/202108014_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()

df_samples = df_samples[grepl("P-norm|D-norm|T-norm|D-COPD", df_samples$condition),]
df_samples = df_samples[!grepl("UNC_44_P|VU_29_D|VU_35_D", df_samples$sample),]


#============ B - Cell numbers in each sample ============================
meta.data = object@meta.data
Cell_types <- c("Cell_subtype","Cell_type","UMAP_land","Family","Superfamily")
cell_number = lapply(Cell_types, function(Cell_type){
    as.data.frame.matrix(table(meta.data[,Cell_type], meta.data$orig.ident))
}) %>% bind_rows()
# remove En UMAP_land
cell_number[grep("En\\...",rownames(cell_number)),]
cell_number = cell_number[-grep("En\\...",rownames(cell_number))[1],]
Symbol = gsub("\\..*","",rownames(cell_number))
cell_number = cell_number[!duplicated(Symbol),]
rownames(cell_number) %<>% gsub("\\..*","",.)

#====
df_number <- readxl::read_excel("doc/Table I. Cell ontology and distribution - template.xlsx",
                                    sheet = "B - Cell numbers in each sample")
colnames(df_number) = df_number[2,]
df_number = df_number[-2,]
df_number = df_number[!is.na(df_number$Symbol),]
Symbol = df_number$Symbol = gsub(" .*","",df_number$Symbol)
wix_columns = colnames(df_number)[-c(1:5)]
setdiff(rownames(cell_number), Symbol) #"TNK"   "BC-S"  "SM-Pr"
setdiff(Symbol,rownames(cell_number)) #

orig.ident =  colnames(cell_number)
rownames(df_number) = Symbol
table(colnames(df_number)[6:ncol(df_number)] == colnames(cell_number))

write.xlsx(cell_number[Symbol,], file = paste0("output/20211007/","B_Cell_numbers.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

meta.data %<>% dplyr::filter(Cell_subtype != "Un")
total_cell_number = as.data.frame(table(meta.data$orig.ident))
table(total_cell_number$Var1 == colnames(cell_number))
colnames(total_cell_number) = c("patient","total")
write.xlsx(total_cell_number, file = paste0(path,"0_patient_cell_number.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
#============ A - Cell proportions per group ============================
df_number <- readxl::read_excel("doc/Table I. Cell ontology and distribution - template.xlsx",
                                sheet = "A - Cell proportions per group")
colnames(df_number) = df_number[2,]
df_number = df_number[-2,]
df_number = df_number[!is.na(df_number$Symbol),!is.na(colnames(df_number))]
Symbol = df_number$Symbol = gsub(" .*","",df_number$Symbol)

cell_number = readxl::read_excel(path = paste0("output/20211007/","B_Cell_numbers.xlsx"),)
colnames(cell_number)[1] = "cell_types"
cell_number %<>% pivot_longer(!cell_types,names_to = "orig.ident", values_to = "num")
cell_number$regions = plyr::mapvalues(cell_number$orig.ident,
                                      from = df_samples$sample,
                                      to = df_samples$regions)
cell_number = cell_number[,c("cell_types","orig.ident","regions","num")]
regions = unique(cell_number$regions)
# mean_res % among all cells (per sample)		
mean_res = cell_number %>% 
    group_by(regions,cell_types) %>%
    summarize(num = sum(num,na.rm = TRUE))
mean_res_all1 <- mean_res

(superfamily <- na.omit(unique(df_number$superfamily)))
mean_res_all2 = cell_number %>% dplyr::filter(cell_types %in% superfamily) %>%
    group_by(regions) %>% 
    summarize(total = sum(num,na.rm = TRUE))
mean_res_all1$total = plyr::mapvalues(mean_res_all1$regions,from = mean_res_all2$regions, 
                                 to = mean_res_all2$total) %>% as.numeric()
mean_res_all1$mean = mean_res_all1$num / mean_res_all1$total *100
mean_res_all = mean_res_all1 %>% pivot_wider(!c("num","total"), names_from = "regions", values_from = "mean")
mean_res_all %<>% column_to_rownames(var="cell_types")

# mean_res % within superfamily (per sample)
mean_res_super1 <- mean_res
mean_res_super2 = mean_res_super1 %>% dplyr::filter(cell_types %in% superfamily)
mean_res_super1$total = 0
for(region in regions){
    for(super in superfamily){
        cells = pull(df_number[df_number$superfamily %in% super,"Symbol"])
        for(cell in cells){
            mean_res_super1[mean_res_super1$regions %in% region & mean_res_super1$cell_types %in% cell,
                       "total"] = mean_res_super2[mean_res_super2$regions %in% region &
                                                      mean_res_super2$cell_types %in% super,"num"]
        }
    }
}
mean_res_super1$mean = mean_res_super1$num / mean_res_super1$total *100
mean_res_super = mean_res_super1 %>% pivot_wider(!c("num","total"), names_from = "regions", values_from = "mean")
mean_res_super %<>% column_to_rownames(var="cell_types")

# mean_res % within family (per sample)
df_number = df_number[!is.na(df_number$Family),]

(Family <- na.omit(unique(df_number$Family)))
mean_res_Family1 <- mean_res
mean_res_Family2 = mean_res_Family1 %>% dplyr::filter(cell_types %in% Family)
mean_res_Family1$total = 0
for(region in regions){
    for(f in Family){
        cells = pull(df_number[df_number$Family %in% f,"Symbol"])
        for(cell in cells){
            mean_res_Family1[mean_res_Family1$regions %in% region & mean_res_Family1$cell_types %in% cell,
                       "total"] = mean_res_Family2[mean_res_Family2$regions %in% region &
                                                  mean_res_Family2$cell_types %in% f,"num"]
        }
    }
}
mean_res_Family1$mean = mean_res_Family1$num / mean_res_Family1$total *100
mean_res_Family = mean_res_Family1 %>% pivot_wider(!c("num","total"), names_from = "regions", values_from = "mean")
mean_res_Family %<>% column_to_rownames(var="cell_types")

table(rownames(mean_res_Family) == rownames(mean_res_all))
table(rownames(mean_res_super) == rownames(mean_res_all))

mean_results = bind_cols(mean_res_all[Symbol,regions],
                         mean_res_super[Symbol,regions],
                         mean_res_Family[Symbol,regions],
                         .name_repair = "universal")
write.xlsx(mean_results, file = paste0(path,"A_Cell_proportions.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#============ A - Cell proportions per group ============================
#========= Differences between groups (2-tailed Mann-Whitney test) =====
df_number <- readxl::read_excel("doc/Table I. Cell ontology and distribution - template.xlsx",
                                sheet = "A - Cell proportions per group")
colnames(df_number) = df_number[2,]
df_number = df_number[-2,]
df_number = df_number[!is.na(df_number$Symbol),!is.na(colnames(df_number))]
Symbol = df_number$Symbol = gsub(" .*","",df_number$Symbol)
#orig.ident = colnames(df_number)[-c(1:5)]

cell_number = readxl::read_excel(path = paste0("output/20211007/","B_Cell_numbers.xlsx"),)
colnames(cell_number)[1] = "cell_types"
cell_number = cell_number[complete.cases(cell_number),]
cell_number %<>% pivot_longer(!cell_types,names_to = "orig.ident", values_to = "num")
cell_number$regions = plyr::mapvalues(cell_number$orig.ident,from = df_samples$sample, to = df_samples$regions)
cell_number = cell_number[,c("cell_types","orig.ident","regions","num")]
regions = unique(cell_number$regions)

# wilcox_res % among all cells (per sample)			
wilcox_res = cell_number %>% group_by(orig.ident,cell_types) %>% 
    summarize(regions = regions,
              num = sum(num,na.rm = TRUE))
wilcox_res_all1 <- wilcox_res

wilcox_res_all2 = cell_number %>% dplyr::filter(cell_types %in% superfamily) %>% 
    group_by(orig.ident) %>% 
    summarize(total = sum(num,na.rm = TRUE))
wilcox_res_all1$total = plyr::mapvalues(wilcox_res_all1$orig.ident,from = wilcox_res_all2$orig.ident, 
                                      to = wilcox_res_all2$total) %>% as.numeric()
wilcox_res_all1$mean = wilcox_res_all1$num / wilcox_res_all1$total *100

compare_list <- list(list("proximal", "distal"),#P vs D
                     list("proximal", "terminal"),#	P vs T
                     list("proximal", c("distal","terminal")),#	P vs D+T
                     list("distal", c("proximal", "terminal")),# D vs P+T
                     list("terminal", "distal"), # T vs D
                     list("terminal", c("proximal", "distal")),#	T vs P+D
                     list("COPD", "distal")#	COPD vs D
                     )
names(compare_list) = c("P vs D","P vs T","P vs D+T","D vs P+T","T vs D","T vs D+P","COPD vs D")
single_wilcox.test <- function(df = wilcox_res_all1, group1 = compare_list[[3]][[1]],
                               group2 = compare_list[[3]][[2]]){
    
    df = df[df$regions %in% unique(c(group1,group2)),]
    df1 = df %>% dplyr::filter(regions %in% group1) %>% pivot_wider(!c("num","total","regions"), 
                                                        names_from = "orig.ident", 
                                                        values_from = "mean") %>% 
        column_to_rownames(var="cell_types")
    df2 = df %>% dplyr::filter(regions %in% group2) %>% pivot_wider(!c("num","total","regions"), 
                                                        names_from = "orig.ident", 
                                                        values_from = "mean") %>%
        column_to_rownames(var="cell_types")
    sdf1 <- split(df1, rownames(df1))
    sdf2 <- split(df2, rownames(df2))
    
    p_value = mapply(function(x,y) {
        scores = wilcox.test(x= unlist(x), y= unlist(y), alternative = "two.sided",exact=FALSE,
                             na.action = na.omit)
        return(scores$p.value)
        }, sdf1, sdf2)
    
    return(p_value)
}

wilcox_res_all = pbapply::pbsapply(compare_list,function(pair) single_wilcox.test(wilcox_res_all1, group1 = pair[[1]], group2 = pair[[2]]))

wilcox_res_all3 = pivot_wider(wilcox_res_all1, !c("regions","num","total"),
                                names_from = orig.ident,values_from = mean)
wilcox_res_all3 %<>% column_to_rownames("cell_types")
write.xlsx(wilcox_res_all3[Symbol,orig.ident], file = paste0(path,"C_Cell_proportions_all.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))



# wilcox_res % within superfamily (per sample)			
(superfamily <- na.omit(unique(df_number$superfamily)))
wilcox_res_super1<- wilcox_res

wilcox_res_super2 = wilcox_res_super1 %>% dplyr::filter(cell_types %in% superfamily)
wilcox_res_super1$total = 0
for(region in regions){
    Samples = df_samples %>% dplyr::filter(regions %in% region) %>%.[,"sample"]
    for(Sample in Samples){
        for(super in superfamily){
            cells = pull(df_number[df_number$superfamily %in% super,"Symbol"])
            for(cell in cells){
                wilcox_res_super1[wilcox_res_super1$orig.ident %in% Sample & wilcox_res_super1$cell_types %in% cell,
                                "total"] = wilcox_res_super2[wilcox_res_super2$orig.ident == Sample &
                                                            wilcox_res_super2$cell_types == super,"num"]
            }
        }
    }
    Progress(which(regions %in% region), length(regions))
}



wilcox_res_super1$mean = wilcox_res_super1$num / wilcox_res_super1$total *100
wilcox_res_super = pbapply::pbsapply(compare_list,function(pair) single_wilcox.test(wilcox_res_super1, group1 = pair[[1]], group2 = pair[[2]]))
wilcox_res_super3 = pivot_wider(wilcox_res_super1, !c("regions","num","total"),
                                names_from = orig.ident,values_from = mean)
wilcox_res_super3 %<>% column_to_rownames("cell_types")
write.xlsx(wilcox_res_super3[Symbol,orig.ident], file = paste0(path,"D_Cell_proportions.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))



# wilcox_res % within family (per sample)			
df_number = df_number[!is.na(df_number$Family),]
(Family <- na.omit(unique(df_number$Family)))

wilcox_res_Family1<- wilcox_res
wilcox_res_Family2 = wilcox_res_Family1 %>% dplyr::filter(cell_types %in% Family)
wilcox_res_Family1$total = 0

#======
for(region in regions){
    Samples = df_samples %>% dplyr::filter(regions %in% region) %>%.[,"sample"]
    for(Sample in Samples){
        for(f in Family){
            cells = pull(df_number[df_number$Family %in% f,"Symbol"])
            for(cell in cells){
                wilcox_res_Family1[wilcox_res_Family1$orig.ident == Sample &
                                       wilcox_res_Family1$cell_types == cell,
                                  "total"] = wilcox_res_Family2[wilcox_res_Family2$orig.ident == Sample &
                                                                    wilcox_res_Family2$cell_types == f,"num"]
            }
        }
    }
    Progress(which(regions %in% region), length(regions))
}

wilcox_res_Family1$mean = wilcox_res_Family1$num / wilcox_res_Family1$total *100
wilcox_res_Family = pbapply::pbsapply(compare_list,function(pair) single_wilcox.test(wilcox_res_Family1, group1 = pair[[1]], group2 = pair[[2]]))
wilcox_res_Family3 = pivot_wider(wilcox_res_Family1, !c("regions","num","total"),
                                names_from = orig.ident,values_from = mean)
wilcox_res_Family3 %<>% column_to_rownames("cell_types")
write.xlsx(wilcox_res_Family3[Symbol,orig.ident], file = paste0(path,"E_Cell_proportions.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))


table(rownames(wilcox_res_super) == rownames(wilcox_res_all))
table(rownames(wilcox_res_Family) == rownames(wilcox_res_all))
table(Symbol %in% rownames(wilcox_res_all))
Symbol[!Symbol %in% rownames(wilcox_res_all)]
wilcox_results = bind_cols(as.data.frame(wilcox_res_all[Symbol,names(compare_list)]),
                           as.data.frame(wilcox_res_super[Symbol,names(compare_list)]),
                           as.data.frame(wilcox_res_Family[Symbol,names(compare_list)]),
                           .name_repair = "universal")
dim(wilcox_results)
write.xlsx(wilcox_results, file = paste0(path,"A_Cell_proportions_wilcox_res.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

