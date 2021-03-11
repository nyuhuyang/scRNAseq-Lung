invisible(lapply(c("Seurat","dplyr","cowplot","kableExtra",
                   "magrittr","eulerr","biomaRt"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

Single_cell <- readRDS(file = "data/Lung_30_20200710.rds") 
load("data/Lung_bulk_20210226.Rda")
bulk <- object
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
