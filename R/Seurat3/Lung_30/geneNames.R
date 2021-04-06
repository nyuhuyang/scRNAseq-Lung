invisible(lapply(c("Seurat","dplyr","cowplot","kableExtra",
                   "magrittr","eulerr","biomaRt"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
set.seed(101)
"Notation:"
"grch37, grch38: gene names from Ensembl datasets"
"hg19, hg38: gene names from single cell or GTEx, with unknown manipulation"

#=========== load grch37 genes ==========
grch37 = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",GRCh = 37)
attributes <-listAttributes(grch37)
head(attributes[,c("name","description")],5)
symbols_grch37 <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), 
                 #filters = 'ensembl_gene_id', 
                 #values = hg_38$ensembl_gene_id, 
                 mart = grch37)
colnames(symbols_grch37)[2] = "grch37"

#=========== read 10X genes ============
list.dirs("data")

sample_id <- c("UNC_52_D","UNC_66_D","VU_37-D")
tsv_list <- pbapply::pblapply(sample_id, function(x){
        file_name = paste0("data/",x,"/outs/filtered_feature_bc_matrix/features.tsv")
        tmp <- read.table(file = gzfile(paste0(file_name,".gz")),
                               sep = '\t', header = F)
        tmp
})
hg_19 <- read.table(file = 'data/UNC_44_D/outs/filtered_gene_bc_matrices/hg19/genes.tsv',
                           sep = '\t', header = F)
colnames(hg_19) = c("ensembl_gene_id","hg19")
hg_19$hg19 %<>% make.unique()
hg_19$hg19 %<>% gsub("_","-",.)
hg_19$single_cell = "single_cell"
# load single-cell
Single_cell <- readRDS(file = "data/Lung_30_20200710.rds") 
single_cell_genes <- rownames(Single_cell) 

table(single_cell_genes %in% symbols_grch37$grch37)
table(single_cell_genes %in% hg_19$hg19)

single_cell_genes[!single_cell_genes %in% symbols_grch37$grch37]

"grch37 doesn't have all gene names from 10X genes,
because make.unique fabricated gene names"

hg_19_37 <- full_join(symbols_grch37, hg_19,by="ensembl_gene_id")
hg_19_37[is.na(hg_19_37)] = ""
missing_37 <- hg_19_37$grch37 %in% hg_19_37[duplicated(hg_19_37$grch37),"grch37"] & 
        stringr::str_length(hg_19_37$hg19) > 0 &
        hg_19_37$grch37 != hg_19_37$hg19
table(missing_37)
hg_19_37[missing_37,"grch37"] = hg_19_37[missing_37,"hg19"]

table(single_cell_genes %in% hg_19_37$hg19)
table(single_cell_genes %in% hg_19_37$grch37)
still_missing <- hg_19_37$hg19 %in% single_cell_genes[!single_cell_genes %in% hg_19_37$grch37]
hg_19_37[hg_19_37$hg19 %in% single_cell_genes[!single_cell_genes %in% hg_19_37$grch37],]

#counts = read.csv("data/RNA-seq/GTEx-Lung-counts.csv",row.names = 1)
tpm = read.csv("data/RNA-seq/GTEx-Lung-tpm.csv",row.names = 1)
# change gene name
hg_38 = data.frame(ensembl_gene_id = gsub("\\.\\d+$","",rownames(tpm)),
                   GTEx = "GTEx",
                   row.names = rownames(tpm))

grch38 = useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <-listAttributes(grch38)
head(attributes[,c("name","description")],20)
symbols_grch38 <- getBM(attributes=c(attributes$name[1],'external_gene_name'), 
                 #filters = 'ensembl_gene_id', 
                 #values = hg_38$ensembl_gene_id, 
                 mart = grch38)
hg_38 %<>% full_join(symbols_grch38,by="ensembl_gene_id")


colnames(hg_38)[3] = "grch38"
hg_19_38 <- full_join(hg_19_37, hg_38,by="ensembl_gene_id")
dim(hg_19_38)
hg_19_38[is.na(hg_19_38)] = ""
hg_19_38 %<>% filter(grch37 != "" | hg19 != "" | grch38 != "")
dim(hg_19_38) # 76292
saveRDS(hg_19_38,"data/RNA-seq/hg_19_38.rds")

# load bulk and replace names
load("data/Lung_bulk_20210226.Rda")
keep_genes <- (rownames(object)  %in% hg_19_38$grch37)
table(keep_genes)
bulk_genes <- rownames(object)[keep_genes]
object %<>% subset(features = bulk_genes)

hg_19_38 = readRDS("data/RNA-seq/hg_19_38.rds")

hg_19_38_short <- hg_19_38 %>% filter(hg19 %in% bulk_genes) %>%
        filter(grch38 != "") %>%
        filter(hg19 != grch38)
dim(hg_19_38_short)
newnames <- plyr::mapvalues(bulk_genes,
                            from = hg_19_38_short$hg19,
                            to = hg_19_38_short$grch38)
table(duplicated(newnames))
newnames %<>% make.unique()
object %<>% RenameGenesSeurat(newnames = newnames)
object %<>% FindVariableFeatures
object@assays$RNA@scale.data = matrix(0,0,0)

rownames(object)[!(rownames(object) %in% c(hg_19_38$grch38,hg_19_38$hg19))]
saveRDS(object, file = "data/Lung_bulk_20210402.rds")
bulk <- readRDS(file = "data/Lung_bulk_20210402.rds")
# replace names for single cell
hg_19_38 = readRDS("data/RNA-seq/hg_19_38.rds")

hg_19_38_short <- hg_19_38 %>% filter(hg19 %in% single_cell_genes) %>%
        filter(grch38 != "") %>%
        filter(hg19 != grch38)
dim(hg_19_38_short)
newnames <- plyr::mapvalues(single_cell_genes,
                            from = hg_19_38_short$hg19,
                            to = hg_19_38_short$grch38)
newnames %<>% make.unique()
Single_cell %<>% RenameGenesSeurat(newnames = newnames)
Single_cell %<>% FindVariableFeatures
Single_cell@assays$SCT@scale.data = matrix(0,0,0)
Single_cell@meta.data$cell_types %<>% gsub("^d-S$","TASC",.)
table(rownames(Single_cell) %in% hg_19_38$grch38)
table(rownames(Single_cell) %in% hg_19_38$hg19)

rownames(Single_cell)[!(rownames(Single_cell) %in% hg_19_38$hg19)]
rownames(Single_cell)[!(rownames(Single_cell) %in% c(hg_19_38$grch38,hg_19_38$hg19))]

Single_cell@assays$RNA = NULL
Single_cell@assays$SCT@scale.data = matrix(0,0,0)
format(object.size(Single_cell), unit = "GB")

saveRDS(Single_cell, file = "data/Lung_SCT_30_20200710.rds")






object <- readRDS(file = "data/Lung_SCT_30_20200710.rds")
GTEx <- data.table::fread("data/RNA-seq/GTEx-Lung-tpm~.csv")
gene_names <- read.csv("output/20210302/gene_name.csv")
Int_supersignatures <- readxl::read_excel("doc/Integrated-signatures-categorized.xlsx") 
signatures = Int_supersignatures[!(Int_supersignatures$Category %in% "HPA"),"Gene"]
signatures = unique(signatures$Gene)
geneNames <- list("Single_cell" = rownames(Single_cell),
                   "bulk" = rownames(bulk),
                   "GTEx" = GTEx$V1,
                  "V1"= unique(gene_names$Single.cell),
                  "V2"= unique(gene_names$GTEx.v8),
                  "signatures"= signatures)
plot(venn(geneNames))
plot(euler(geneNames, shape = "circle"), quantities = TRUE)
plot(euler(geneNames[c("Single_cell","GTEx","signatures")], shape = "circle"), quantities = TRUE)
plot(euler(geneNames[c("Single_cell","GTEx","signatures","V1")], shape = "circle"), quantities = TRUE)
plot(euler(geneNames[c("Single_cell","signatures","V1")], shape = "circle"), quantities = TRUE)
plot(euler(geneNames[c("GTEx","V2")], shape = "circle"), quantities = TRUE)

replace_gene_names <- read.csv("output/20210302/replace_gene_names.csv")
Int_supersignatures <- readxl::read_excel("doc/Integrated-signatures-categorized.xlsx") 

Int_supersignatures$Gene %<>% plyr::mapvalues(from = replace_gene_names$Single.cell,
                                              to = replace_gene_names$GTEx.v8)
openxlsx::write.xlsx(Int_supersignatures, file =  "doc/Integrated-signatures-categorized-newNames.xlsx",
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#==================
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
datasets[datasets$dataset %in% "hsapiens_gene_ensembl",]
GRCh38 <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
head(listAttributes(GRCh38))
GRCh38_genes <- getBM(attributes=c("ensembl_gene_id",'hgnc_symbol'), 
                      mart = GRCh38)


listEnsemblArchives()
ensembl_archive <- useMart("ensembl", host="http://dec2016.archive.ensembl.org/grch37/")
ensembl_archive <- useEnsembl(biomart = "ensembl", 
                              dataset = "hsapiens_gene_ensembl", 
                              version = "87",
                              GRCh = "37")
listDatasets(ensembl_archive)
GRCh37 <- useDataset("hsapiens_gene_ensembl",mart=ensembl_archive)
head(listAttributes(GRCh37))
GRCh37_genes <- getBM(attributes=c("ensembl_gene_id",'hgnc_symbol'), 
                      mart = GRCh37)

geneNames <- list("Single_cell" = rownames(Single_cell),
                  "bulk" = rownames(bulk),
                  "GTEx" = GTEx$V1,
                  "V1"= unique(gene_names$Single.cell),
                  "V2"= unique(gene_names$GTEx.v8),
                  "signatures"= signatures,
                  "GRCh37" = unique(GRCh37_genes$hgnc_symbol),
                  "GRCh38" = unique(GRCh38_genes$hgnc_symbol))
plot(euler(geneNames[c("GTEx","GRCh38")], shape = "circle"), quantities = TRUE)
plot(euler(geneNames[c("GRCh37","Single_cell")], shape = "circle"), quantities = TRUE)

plot(euler(geneNames[c("bulk","GRCh38")], shape = "circle"), quantities = TRUE)
plot(euler(geneNames[c("Single_cell","GRCh38","GRCh37")], shape = "circle"), quantities = TRUE)



# GTEX
genes <-scan("../seurat_resources/Enrichr/GTEx_Tissue_Sample_Gene_Expression_Profiles_up.txt", what="", sep="\t")
genes = genes[-grep(" ",genes)] %>% unique
table(genes %in% hg_19_38$hg19)
table(genes %in% hg_19_38$hg38)

genes[genes %in% hg_19_38$hg19]
genes[!(genes %in% hg_19_38$hg38)]
