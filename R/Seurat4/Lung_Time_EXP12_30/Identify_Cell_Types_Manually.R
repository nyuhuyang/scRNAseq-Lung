library(Seurat)
library(dplyr)
library(tidyr)
#library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
library(stringr)
library(openxlsx)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load data
meta.data <- readRDS("output/Lung_time15_metadata_20220523.rds")
meta.data$sample <- meta.data$patient
#======== rename ident =================
df_annotation <- readxl::read_excel("doc/Annotations/20230125 Annotation ALI EX-12-30 RS.xlsx",
                                    sheet = "annotation")
df_annotation$output <- "unknown"
resolutions = c("SCT_snn_res.0.8","SCT_snn_res.2")
meta.data$`cell state` = "unknown"
for(i in 1:length(resolutions)){
    resolution_ind = which(!is.na(pull(df_annotation[,resolutions[i]])))
        for(m in resolution_ind){
            cl = pull(df_annotation[m,resolutions[i]])
            change_to <- pull(df_annotation[m,"cell state"])
            select_id <- meta.data %>% dplyr::filter(!!as.name(resolutions[i]) %in% cl)
            if(!is.na(df_annotation$`subsetting step`[m])) {
                subsetting <- df_annotation$`subsetting step`[m]
                select_id <- dplyr::filter(select_id,eval(parse(text = subsetting)))
            } else subsetting <- ""

            meta.data[rownames(select_id),"cell state"] <- change_to
            output <- paste(nrow(select_id),"cells in",resolutions[i],"at cluster",cl,subsetting,"------->",change_to)
            print(output)
            df_annotation[m,"output"] <- output
        }
}
openxlsx::write.xlsx(df_annotation, file =  paste0(path,"20230125 Annotation ALI EX-12 30 YH.xlsx"),
                     colNames = TRUE,rowNames = FALSE,borders = "surrounding")

Cell_types <- c("cell state","cell type","cell group","Stage","Path")

df_annotation %<>% as.data.frame()
df_annotation <- df_annotation[!duplicated(df_annotation$`cell state`),]
for(Cell_type in Cell_types){
    meta.data[,Cell_type] = plyr::mapvalues(meta.data$`cell state`,
                                            from = df_annotation[,"cell state"],
                                            to = df_annotation[,Cell_type])
}

saveRDS(meta.data, file = "output/Lung_time15_metadata_20220523_v2.rds")
