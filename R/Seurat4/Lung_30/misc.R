# conda activate r4.0.3
library(Seurat)
library(data.table)
library("openxlsx")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
object %<>% subset(subset = Family == "En" &
                       Regions %in% c("COPD","distal"))
df = FetchData(object, vars = c("F2R","Cell_type","Regions"))
df_list <- split(df,f = df$Cell_type)
df_list[["En"]] = df
df_list[["En-SM"]] = NULL
df_list = df_list[c("En","En-v","En-a","En-c","En-l")]
write.xlsx(df_list, file = paste0(path,"F2R_expression_En.xlsx"),
           colNames = TRUE, rowNames = TRUE,borders = "surrounding")


meta.data = readRDS("~/Downloads/annot_umap.rds")
meta.data = meta.data[,c("cluster_id","sub_cluster")]
exp = readRDS("~/Downloads/expression_values.rds")
genes = intersect(rownames(object),colnames(exp))
length(genes)
exp = as.data.frame(exp[,genes])
df = data.table(cbind(meta.data, exp))
df = df[order(df$sub_cluster),]
fwrite(df,file = paste0(path,"human_embryo.csv"), row.names=FALSE, quote=FALSE)
test = t(df)
fwrite(test,file = paste0(path,"human_embryo_trans.csv"), row.names=FALSE,
       col.names = FALSE, quote=FALSE)

test = cbind("gene" = rownames(test),test)
write.xlsx(df, file = paste0(path,"human_embryo.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

genes1 = intersect(VariableFeatures(object),colnames(exp))
length(genes1)
exp1 = as.data.frame(exp[,genes1])
df1 = cbind(meta.data, exp1)
df1 = df1[order(df1$sub_cluster),]

fwrite(df1,file = paste0(path,"human_embryo_2000_genes.csv"), row.names=FALSE, quote=FALSE)

write.xlsx(df1, file = paste0(path,"human_embryo_2000_genes.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
