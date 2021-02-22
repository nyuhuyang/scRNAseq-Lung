library(CePa)
library(biomaRt)
library(dplyr)
library(magrittr)
library(Seurat)
library(kableExtra)
library(openxlsx)
library(data.table)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
set.seed(101)
# on linux
GTE_meta.data = read.csv("data/RNA-seq/GTEx Portal - lung sample info.csv",
                         stringsAsFactors = F)
rownames(GTE_meta.data) = GTE_meta.data$Tissue.Sample.ID
rownames(GTE_meta.data) = gsub("-","\\.", rownames(GTE_meta.data))
select = paste(rownames(GTE_meta.data),collapse = "|")

# counts
counts = CePa::read.gct("data/RNA-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
format(object.size(counts),unit = "GB")

select = grep(select,colnames(counts), value = T)

counts = counts[,select]
counts = counts[rowSums(counts) > 0,]
write.csv(counts, file="data/RNA-seq/GTEx-Lung-counts.csv")

# tpm
tpm =  CePa::read.gct("data/RNA-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
format(object.size(tpm),unit = "GB")
select = grep(select,colnames(tpm), value = T)

tpm = tpm[,select]
tpm = tpm[rowSums(tpm) > 0,]
write.csv(tpm, file="data/RNA-seq/GTEx-Lung-tpm.csv")
#fwrite(tpm, file="data/RNA-seq/GTEx-Lung-tpm.csv")

# on Mac
GTE_meta.data = read.csv("data/RNA-seq/GTEx Portal - lung sample info.csv",
                         stringsAsFactors = F)
rownames(GTE_meta.data) = GTE_meta.data$Tissue.Sample.ID

#counts = read.csv("data/RNA-seq/GTEx-Lung-counts.csv",row.names = 1)
counts = read.csv("data/RNA-seq/GTEx-Lung-tpm.csv",row.names = 1)
df_counts = data.frame(ensembl_gene_id = gsub("\\..*","",rownames(counts)),
                       row.names = rownames(counts))

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
head(attributes,20)
symbols <- getBM(attributes=c("ensembl_gene_id",'hgnc_symbol'), 
                filters = 'ensembl_gene_id', 
                values = df_counts$ensembl_gene_id, 
                mart = ensembl)
df_counts$gene <- symbols[match(df_counts$ensembl_gene_id, symbols$ensembl_gene_id),
                      "hgnc_symbol"]
df_counts %<>% cbind(counts)
df_counts = subset(df_counts, gene != "")
df_counts = df_counts[!duplicated(df_counts$gene),]
rownames(df_counts) = df_counts$gene
df_counts = df_counts[, -grep("ensembl_gene_id|gene",colnames(df_counts))]
df_counts %<>% as.matrix()

colnames(df_counts) %<>% gsub("\\.SM.*","",.)
rownames(GTE_meta.data) %<>% gsub("-",".",.)
GTE_meta.data = GTE_meta.data[colnames(df_counts),]
table(rownames(GTE_meta.data) == colnames(df_counts) )
object <- CreateSeuratObject(counts = df_counts, meta.data = GTE_meta.data,
                             min.cells = 3, min.features = 200)
object %<>% SCTransform
object %<>% FindVariableFeatures(selection.method = "vst",
                               num.bin = 20,nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), 
                               dispersion.cutoff = c(1, Inf))
object %<>%  RunPCA
object %<>%  FindNeighbors(dims = 1:20)
object %<>%  FindClusters
object %<>%  RunUMAP(dims = 1:20)
object@meta.data$Age = gsub("-.*","",object$Age.Bracket) %>% as.integer()
object@meta.data$SCT_snn_res.0.8 %<>% as.numeric()
object@meta.data %>%
        group_by(SCT_snn_res.0.8, Age, Sex) %>%
        summarize(n())
object@meta.data %>%
        group_by(SCT_snn_res.0.8) %>%
        summarize(mean = mean(Age))
UMAPPlot.1(object, group.by = "Age.Bracket", do.print = T, cols = Singler.colors)#Sex Age.Bracket
FeaturePlot.1(object,features = "nCount_RNA")
FeaturePlot.1(object,features = "Age")
save(object, file = "data/Lung_GTEx_20200307.Rda")

# write counts
df_counts  = as.matrix(object[["RNA"]]@counts)
df_counts = cbind(object[["Age.Bracket"]],t(df_counts))
df_counts = df_counts[order(df_counts[,"Age.Bracket"]),]
write.csv(t(df_counts), file = "data/RNA-seq/GTEx-Lung-tpm.csv")
df_counts = read.csv("data/RNA-seq/GTEx-Lung-counts.csv",row.names = 1)

# write TPM with  final integrated signatures
library(tibble)
tpm = read.csv("data/RNA-seq/GTEx-Lung-tpm.csv",row.names = 1)
tpm = tpm[-1,]
tpm %<>% rownames_to_column(var = "gene")
tpm[,2:ncol(tpm)] %<>% apply(2,as.numeric)
df = readxl::read_excel("doc/Integrated-signatures-categorized.xlsx")
colnames(df) %<>% tolower()
sheet = unique(tpm$`excel sheet`)
tpm %<>% left_join(df,., by = "gene")
tpm_list <- split(tpm,f = tpm$`excel sheet`)
tpm_list = tpm_list[sheet]
tpm_list %<>% lapply(function(x) x[,-4])
write.xlsx(tpm_list, file = "Yang/GTEx/Integrated-signatures-categorized - 577 GTEx samples ordered.xlsx",
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))



# read optimized DEGs
DEGs <- read.csv(file = "Yang/Lung_30/DE_analysis/Optimized_cell_type_DEG.csv",row.names = 1)
DEGs = DEGs[,c("gene","cluster")]
colnames(DEGs)[2] ="cell.type"

#  ======== DE Between different age groups============
(load(file = "data/Lung_GTEx_20200307.Rda"))
Idents(object) = "Age.Bracket"
object %<>% sortIdent()
gender_marker <- FindAllMarkers.UMI(object, p.adjust.methods = "BH",
                                    logfc.threshold = 0,
                                    return.thresh = 1,
                                    only.pos = F, 
                                    min.pct = 0.1)
write.csv(gender_marker,paste0(path,"age_markers_FC0.csv"))


#  ======== DE Between different age groups, separately============
(load(file = "data/Lung_GTEx_20200307.Rda"))
Young_ages <- list(c("20-29","30-39"),
                   c("20-29","30-39"),
                   c("20-29","30-39"),
                   c("20-29"),
                   c("20-29"),
                   c("20-29"),
                   c("20-29"))
old_ages <- list(c("60-69","70-79"),
                 c("50-59","60-69","70-79"),
                 c("40-49","50-59","60-69","70-79"),
                 c("60-69","70-79"),
                 c("50-59","60-69","70-79"),
                 c("40-49","50-59","60-69","70-79"),
                 c("30-39","40-49","50-59","60-69","70-79"))

Idents(object) = "Age.Bracket"
object %<>% sortIdent()
table(Idents(object))
age_markers <- FindPairMarkers(object, 
                               ident.1 = Young_ages,
                               ident.2 = old_ages,
                               p.adjust.methods = "BH",
                               logfc.threshold = 0,
                               return.thresh = 1,
                               only.pos = F, 
                               min.pct = 0.1)

age_markers$FC = 2^(age_markers$avg_logFC)
#age_markers = age_markers[age_markers$p_val_adj<0.05,]
write.csv(age_markers,paste0(path,"age_markers_subset_FC0.csv"))


write.xlsx(age_markers_list, file = paste0(path,"2a-age_gender_markers_FC0.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#  ======== DE Between different age groups, males and females separately============
(load(file = "data/Lung_GTEx_20200307.Rda"))
Idents(object) = "Sex"
Young_ages <- list(c("20-29","30-39"),
                   c("20-29","30-39"),
                   c("20-29","30-39"),
                   c("20-29"),
                   c("20-29"),
                   c("20-29"),
                   c("20-29"))
old_ages <- list(c("60-69","70-79"),
                 c("50-59","60-69","70-79"),
                 c("40-49","50-59","60-69","70-79"),
                 c("60-69","70-79"),
                 c("50-59","60-69","70-79"),
                 c("40-49","50-59","60-69","70-79"),
                 c("30-39","40-49","50-59","60-69","70-79"))
age_markers_list <- list()
for(gender in c("male","female")){
        sub_object = subset(object, idents = gender)
        Idents(sub_object) = "Age.Bracket"
        sub_object %<>% sortIdent()
        table(Idents(sub_object))
        age_markers <- FindPairMarkers(sub_object, 
                                       ident.1 = Young_ages,
                                       ident.2 = old_ages,
                                       p.adjust.methods = "BH",
                                       logfc.threshold = 0,
                                       return.thresh = 1,
                                       only.pos = F, 
                                       min.pct = 0.1)
        
        age_markers$FC = 2^(age_markers$avg_logFC)
        age_markers_list[[gender]] = age_markers
        #age_markers = age_markers[age_markers$p_val_adj<0.05,]
        write.csv(age_markers,paste0(path,"age_markers_",gender,"_FC0.csv"))
}

write.xlsx(age_markers_list, file = paste0(path,"2a-age_gender_markers_FC0.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#  ======== DE Between different age groups, males and females separately============
(load(file = "data/Lung_GTEx_20200307.Rda"))
object$Age_Sex = paste0(object$Age.Bracket, "_", object$Sex)
ages <- list(c("20-29","30-39"),
             c("40-49","50-59","60-69","70-79"),
             c("50-59","60-69","70-79"),
             c("60-69","70-79"))
ages_male = lapply(ages, function(x) paste0(x, "_","male"))
ages_female = lapply(ages, function(x) paste0(x, "_","female"))
Idents(object) = "Age_Sex"
age_markers <- FindPairMarkers(object, 
                               ident.1 = ages_male,
                               ident.2 = ages_female,
                               p.adjust.methods = "BH",
                               logfc.threshold = 0,
                               return.thresh = 1,
                               only.pos = F, 
                               min.pct = 0.1)
age_markers$FC = 2^(age_markers$avg_logFC)
#age_markers = age_markers[age_markers$p_val_adj<0.05,]
write.csv(age_markers,paste0(path,"2b-age_gender_markers_FC0.csv"))

Idents(object) = "Sex"
gender_marker <- FindAllMarkers.UMI(object, p.adjust.methods = "BH",
                                    logfc.threshold = 0,
                                    return.thresh = 1,
                                    only.pos = F, 
                                    min.pct = 0.1)
write.csv(gender_marker,paste0(path,"2b-gender_markers_FC0.csv"))

#=======

age_markers <- read.csv(paste0("Yang/GTEx/","age_markers_FC0.csv"),row.names = 1)

age_markers %>% group_by()
(top <-  age_markers %>% 
        group_by(cluster1.vs.cluster2) %>% 
        top_n(2, avg_logFC))
grep("SCGB3A2", age_markers$gene)


GCT <- c("Ion-NEC-morpheus-output-re","C1-C4-H-morpheus-output-re.gct")
for(i in seq_along(GCT)){
        counts = read.gct(paste0("data/RNA-seq/",GCT[i],".gct"))
        colnames(counts) = gsub("\\.","-",colnames(counts))
        write.csv(counts,paste0(path,GCT[i],".csv"))
}
#  ======== expression ============
(load(file = "data/Lung_GTEx_20200307.Rda"))
meta.data = object@meta.data
meta.data = meta.data[,c("Sex","Age.Bracket" ,"Hardy.Scale","Pathology.Notes")]

meta.data = meta.data[order(meta.data$Age.Bracket, decreasing = F),]
meta.data = meta.data[order(meta.data$Sex, decreasing = T),]
counts <- object@assays$SCT@data %>% as.matrix() %>% t
counts_meta <- cbind(meta.data, counts[rownames(meta.data),])
write.csv(t(counts_meta), paste0(path, "GTEx_counts.csv"))


tpm = read.csv("data/RNA-seq/GTEx-Lung-tpm.csv",row.names = 1)
tpm = tpm[-1,]
tpm %<>% tibble::rownames_to_column(var = "gene")

EVGs.path <- "Yang/Lung_30/DE_analysis/C_Cell_types/"
DEGs = readxl::read_excel(paste0(EVGs.path,"supersignatures.xlsx"))
colnames(DEGs) %<>% gsub("cluster","cell.type",.)
DEGs = DEGs[,c("cell.type","gene")]
DEGs = DEGs[!duplicated(DEGs$gene),]
tpm %<>% right_join(DEGs,by = "gene")
tpm %<>% tibble::column_to_rownames(var = "gene")
tpm = tpm[,-grep("cell.type",colnames(tpm))]
write.csv(tpm, paste0(path, "GTEx_TPM_supersignatures.csv"))

library(dplyr)
library(magrittr)
library(openxlsx)
# assemble a file containing cell type genes with FC and p adj values retrieved from the April analysis. 

# read optimized DEGs
GTEx.path = "Yang/GTEx/results_2020April/"
if(GTEx.path == "Yang/GTEx/results_2020April/"){
        DEGs <- read.csv(file = "Yang/Lung_30/DE_analysis/Optimized_cell_type_DEG.csv",row.names = 1)
        DEGs = DEGs[,c("gene","cluster")]
        colnames(DEGs)[2] ="cell.type"
}

# read optimized DEGs
GTEx.path = "Yang/GTEx/results_2020December/"
if(GTEx.path == "Yang/GTEx/results_2020December/"){
        # Age-related GTEx analysis – you would need to extract only data for EVG genes from age-related GTEx genes
        EVGs.path <- "Yang/Lung_30/DE_analysis/F_EVGs_allCells/"
        superfamily <- c("Epithelial","Structural","Immune")
        EVGs_list <- lapply(superfamily, function(s) {
                tmp = readxl::read_excel(paste0(EVGs.path,"Lung_30-EVGs.xlsx"), sheet = s)
                tmp[,-1]
        })
        cbind.fill <- function(...){
                nm <- list(...) 
                nm <- lapply(nm, as.matrix)
                n <- max(sapply(nm, nrow)) 
                do.call(cbind, lapply(nm, function (x) 
                        rbind(x, matrix("", n-nrow(x), ncol(x))))) 
        }
        EVGs <- do.call(cbind.fill,EVGs_list)
        EVGs <- tidyr::gather(as.data.frame(EVGs),cell.type, gene)
        EVGs = EVGs[!is.na(EVGs$gene),]
        DEGs = EVGs
}
# 20210209 =============
# read supersignatures
GTEx.path = "Yang/GTEx/results_2021Feb/"
if(!dir.exists(GTEx.path))dir.create(GTEx.path, recursive = T)

# Age-related GTEx analysis – you would need to extract only data for EVG genes from age-related GTEx genes
read.path <- "Yang/Lung_30/DE_analysis/C_Cell_types/"
DEGs = read.csv(paste0(read.path,"supersignatures_Extended.csv"))
colnames(DEGs) %<>% gsub("cluster","cell.type",.)
DEGs = DEGs[,c("cell.type","gene")]
# combine with previous results_2020April
read.path = "Yang/GTEx/Archive/"

#1a. age_markers
dge = read.csv(paste0(read.path, "age_markers_FC0.csv"))
if("X" %in% colnames(dge)) dge = dge[, -which(colnames(dge) %in% "X")]
dge %<>% right_join(DEGs, by = "gene") 
dge = dge[order(dge$cell.type),]
rownames(dge) = NULL
write.csv(dge, file = paste0(GTEx.path,"1a.All (both males and females) young vs old.csv"))

#1b. age_markers
dge = read.csv(paste0(read.path, "age_markers_subset_FC0.csv"))
if("X" %in% colnames(dge)) dge = dge[, -which(colnames(dge) %in% "X")]
dge %<>% right_join(DEGs, by = "gene") 
dge = dge[order(dge$cell.type),]
write.csv(dge, file = paste0(GTEx.path,"1b.All (both males and females) part of young vs old.csv"))


#2a. age_gender_markers
for(i in 1:2){
        sex = c("male","female")[i]
        dge <- readxl::read_excel(paste0(read.path,"2a-age_gender_markers_FC0.xlsx"),sheet = sex)
        dge %<>% right_join(DEGs, by = "gene")
        dge = dge[order(dge$cell.type),]
        write.xlsx(dge, file = paste0(GTEx.path,"2a.",sex,"-young_vs_old.xlsx"),
                   colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
        svMisc::progress(i, 2)
}

#"2b-age_gender_markers"
dge <- read.csv(paste0(read.path,"2b-age_gender_markers_FC0.csv"), row.names = 1)
dge %<>% right_join(DEGs, by = "gene")
dge = dge[order(dge$cell.type),]
write.csv(dge, file = paste0(GTEx.path,"2b-male_vs_female_by_age.xlsx"))

dge <- read.csv(paste0(read.path,"2b-gender_markers_FC0.csv"), row.names = 1)
dge %<>% right_join(DEGs, by = "gene")
dge = dge[order(dge$cell.type),]
write.csv(dge, file = paste0(GTEx.path,"2b-male_vs_female.csv"))

# 20210222 =============
# read supersignatures
GTEx.path = "Yang/GTEx/results_2021Feb~/"
if(!dir.exists(GTEx.path))dir.create(GTEx.path, recursive = T)

df = readxl::read_excel("doc/Integrated-signatures-categorized.xlsx")
colnames(df) %<>% tolower()
DEGs = df[,c("gene","family","category")]
# combine with previous results_2020April
read.path = "Yang/GTEx/Archive/"

#1a. age_markers
dge = read.csv(paste0(read.path, "age_markers_FC0.csv"))
if("X" %in% colnames(dge)) dge = dge[, -which(colnames(dge) %in% "X")]
dge %<>% left_join(DEGs, ., by = "gene") 
rownames(dge) = NULL
write.csv(dge, file = paste0(GTEx.path,"1a.All (both males and females) young vs old.csv"))

#1b. age_markers
dge = read.csv(paste0(read.path, "age_markers_subset_FC0.csv"))
if("X" %in% colnames(dge)) dge = dge[, -which(colnames(dge) %in% "X")]
dge %<>% left_join(DEGs,., by = "gene") 
#dge = dge[order(dge$cell.type),]
write.csv(dge, file = paste0(GTEx.path,"1b.All (both males and females) part of young vs old.csv"))


#2a. age_gender_markers
for(i in 1:2){
        sex = c("male","female")[i]
        dge <- readxl::read_excel(paste0(read.path,"2a-age_gender_markers_FC0.xlsx"),sheet = sex)
        dge %<>% left_join(DEGs,., by = "gene")
        #dge = dge[order(dge$cell.type),]
        write.xlsx(dge, file = paste0(GTEx.path,"2a.",sex,"-young_vs_old.xlsx"),
                   colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
        svMisc::progress(i, 2)
}

#"2b-age_gender_markers"
dge <- read.csv(paste0(read.path,"2b-age_gender_markers_FC0.csv"), row.names = 1)
dge %<>% left_join(DEGs,., by = "gene")
#dge = dge[order(dge$cell.type),]
write.csv(dge, file = paste0(GTEx.path,"2b-male_vs_female_by_age.xlsx"))

dge <- read.csv(paste0(read.path,"2b-gender_markers_FC0.csv"), row.names = 1)
dge %<>% left_join(DEGs,., by = "gene")
#dge = dge[order(dge$cell.type),]
write.csv(dge, file = paste0(GTEx.path,"2b-male_vs_female.csv"))
