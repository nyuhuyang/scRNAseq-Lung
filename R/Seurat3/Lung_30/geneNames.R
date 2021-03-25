invisible(lapply(c("Seurat","dplyr","cowplot","kableExtra",
                   "magrittr","eulerr","biomaRt"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
set.seed(101)
#=========== read 10X genes ============
list.dirs("data")

sample_id <- c("UNC_52_D","UNC_66_D","VU_37-D")
tsv_list <- pbapply::pblapply(sample_id, function(x){
        file_name = paste0("data/",x,"/outs/filtered_feature_bc_matrix/features.tsv")
        tmp <- read.table(file = gzfile(paste0(file_name,".gz")),
                               sep = '\t', header = F)
        tmp
})
tsv_list[[4]] <- read.table(file = 'data/UNC_44_D/outs/filtered_gene_bc_matrices/hg19/genes.tsv',
                           sep = '\t', header = F)
identical(tsv_list[[1]],tsv_list[[2]])
identical(tsv_list[[1]],tsv_list[[3]])
tsv_list[[3]]$V3 = NULL
identical(tsv_list[[3]],tsv_list[[4]])
hg_19 <- tsv_list[[4]]
colnames(hg_19) = c("ensembl_gene_id","gene_87")
hg_19$gene_87 %<>% make.unique()
hg_19$gene_87 %<>% gsub("_","-",.)
hg_19$single_cell = "single_cell"
# load single-cell
Single_cell <- readRDS(file = "data/Lung_30_20200710.rds") 
single_cell_genes <- rownames(Single_cell) 
load("data/Lung_bulk_20210226.Rda")
bulk <- object;rm(object);GC()
table(rownames(Single_cell)  %in% hg_19$gene_87)
rownames(bulk)[!(rownames(bulk)  %in% hg_19$gene_87)]

hg_19$gene_87[!(hg_19$gene_87 %in% single_cell_genes)]


#counts = read.csv("data/RNA-seq/GTEx-Lung-counts.csv",row.names = 1)
tpm = read.csv("data/RNA-seq/GTEx-Lung-tpm.csv",row.names = 1)
# change gene name
hg_38 = data.frame(ensembl_gene_id = gsub("\\.\\d+$","",rownames(tpm)),
                   GTEx = "GTEx",
                       row.names = rownames(tpm))

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
head(attributes,20)
symbols <- getBM(attributes=c("ensembl_gene_id",'hgnc_symbol'), 
                 #filters = 'ensembl_gene_id', 
                 #values = hg_38$ensembl_gene_id, 
                 mart = ensembl)
hg_38 %<>% full_join(symbols,by="ensembl_gene_id")
colnames(hg_38)[3] = "gene_103"
hg_19_38 <- full_join(hg_19, hg_38,by="ensembl_gene_id")
dim(hg_19_38)
hg_19_38 %<>% as.matrix()
hg_19_38[hg_19_38[,"gene_103"] == "","gene_103"] = NA
hg_19_38 %<>% as.data.frame
saveRDS(hg_19_38,"data/RNA-seq/hg_19_38.rds")

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
