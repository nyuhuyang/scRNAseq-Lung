path = "Yang/proximal_distal_terminal/Non-Integration/DEGs/cell_types_pairwise/"
SCT_files <- list.files(paste0(path,"SCT"))
RNA_files <- list.files(paste0(path,"RNA"))
SCT_files <- gsub("Lung_24_SCT_pairwise_","",SCT_files)
SCT_files <- gsub("_.*","",SCT_files) %>% as.numeric() %>% sort
num <- 1:80
num[!(num %in% SCT_files)]

RNA_files <- gsub("Lung_24_RNA_pairwise_","",RNA_files)
RNA_files <- gsub("_.*","",RNA_files) %>% as.numeric() %>% sort
num[!(num %in% RNA_files)]

df_labels_Pairs <- readxl::read_excel("doc/DEG cell types comparison plan.xlsx")
Cell.type <- readxl::read_excel("doc/Cell type abbreviation.xlsx")

df_labels_Pairs = df_labels_Pairs[!is.na(df_labels_Pairs$`Group 1`),]
(ident1 <- strsplit(df_labels_Pairs$`Group 1`, "\\+") %>% unlist())
(ident2 <- strsplit(df_labels_Pairs$`Group 2`, "\\+")%>% unlist())
Uni <- c(ident1,ident2)
Uni[!Uni %in% Cell.type$Abbreviation]


# Add Abbreviation
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

conditions = c("proximal","distal","terminal","All")
(con <- conditions[args])

df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")

(load(file = paste0("data/Epi_24_",con,"_20191223.Rda")))
object$cell_types <- plyr::mapvalues(object$cell.types,
                                     from = df_cell_types$`Cell types`,
                                     to = df_cell_types$Abbreviation)
object$cell_types.colors = object$cell.types.colors
Idents(object) = "cell_types"

df <- colData(cds)
if(table(rownames(df) == colnames(object))) {
        df$cell_types = object$cell_types
        df$cell_types.colors = object$cell_types.colors
}
colData(cds) = df
