source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
library(data.table)
library(countToFPKM)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)


counts_df <- read.delim("data/RNA-seq/Counts/SR_EX42_CTR_S1.SR_EX42_IFN_S2.RS58_3_S3.RS58_5_S4.bam.txt",sep = "\t",header = TRUE)
dim(counts_df)
counts_df <- counts_df[complete.cases(counts_df),]
dim(counts_df)
rownames(counts_df) = NULL
counts <- counts_df[,c("Geneid","SR_EX42_CTR_S1","SR_EX42_IFN_S2","RS58_3_S3","RS58_5_S4")] %>% 
    tibble::column_to_rownames("Geneid") %>% as.data.frame %>% as.matrix()

MEAN_INSERT_SIZE = c(262.525491,278.592549,268.5006,295.211964)

fpkm_matrix <- fpkm (counts, counts_df$Length, MEAN_INSERT_SIZE)

write.csv(fpkm_matrix,paste0(path,"new_bulkseq_fpkm.csv"))


rpkm <- function(count, lengths) {
    rate <- count / lengths 
    rate / sum(count) * 1e9 
}

fpkm_list <- apply(counts,2, function(count= x, lengths= counts_df$Length) rpkm(count,lengths))

commen_genes <- intersect(rownames(fpkm_matrix),rownames(fpkm_list))
fpkm_matrix[tail(commen_genes,10),]
fpkm_list[tail(commen_genes,10),]
#====================================
MergedFPKMs <- data.table::fread("data/RNA-seq/MergedFPKMs all samples - 2-24-22 labeled-DONOR INFO.csv") %>% as.data.frame
Class = sapply(MergedFPKMs,class)
(character_Class <- Class[Class == "character"])
for(col in names(character_Class[2:4])) {
    MergedFPKMs[,col] %<>% gsub("[[:space:]]", "", .) %>% as.numeric()
}

tedious_genes <- data.frame(c("1-Dec","BHLHE40"),
                            c("1-Mar","RTL1"),#Mar1
                            c("1-Mar","MARCHF1"),#MARCH1
                            c("1-Sep","SEPTIN1"),
                            c("10-Mar","MARCHF10"),
                            c("10-Sep","SEPTIN1"),
                            c("11-Mar","MARCHF11"),
                            c("11-Sep","SEPTIN11"),
                            c("12-Sep","SEPTIN12"),
                            c("14-Sep","SEPTIN14"),
                            c("15-Sep","SELENOF"),
                            c("2-Mar","PEG10"),#Mar2
                            c("2-Mar","MARCHF2"),#MARCH2
                            c("3-Mar","MARCHF3"),#MARCH3
                            c("3-Sep","SEPTIN3"),
                            c("4-Mar","MARCHF4"),#MARCH4?
                            c("4-Sep","SEPTIN4"),
                            c("5-Mar","MARCHF5"),#MARCH5?
                            c("5-Sep","SEPTIN5"),
                            c("6-Mar","MARCHF6"),#MARCH6?
                            c("6-Sep","SEPTIN6"),
                            c("7-Mar","MARCHF7"),#MARCH7?
                            c("7-Sep","SEPTIN7"), 
                            c("8-Mar","MARCHF8"),#MARCH8?
                            c("8-Sep","SEPTIN8"), 
                            c("9-Mar","MARCHF9"),#MARCH9?
                            c("9-Sep","SEPTIN7")) %>% t
for(i in seq_along(rownames(tedious_genes))) {
    MergedFPKMs$genes %<>% sub(paste0("^",tedious_genes[i,1],"$"), tedious_genes[i,2],.)
}
dim(MergedFPKMs)
MergedFPKMs <- MergedFPKMs[complete.cases(MergedFPKMs),]
dim(MergedFPKMs)
markers <- c("SCGB1A1","SCGB3A2","SFTPB","KRT15","SERPINB3","MUC5AC","SCGB1A1","SCGB3A2","SFTPB","FOXA2")
test <- MergedFPKMs[MergedFPKMs$genes %in% markers,]
rownames(test) = NULL
#https://bioinformatics.stackexchange.com/questions/15573/what-is-the-best-way-to-programmatically-convert-ensembl-ids-from-one-release-to
library(biomaRt)
library(dplyr)

mart37 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  host    = "https://grch37.ensembl.org",
                  path    = "/biomart/martservice",
                  dataset = "hsapiens_gene_ensembl")

mart38 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  #path    = "/biomart/martservice",
                  dataset = "hsapiens_gene_ensembl")
out37 <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id','clone_based_vega_transcript'),
               mart = mart37)
out37 <- out37[!duplicated(out37$ensembl_gene_id),]
colnames(out37) = c("hg19","ensembl_gene_id","hg19_transcript")
out37$hg19_transcript %<>% gsub("-[0-9][[0-9][[0-9]$","",.)

out38 <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id'),
               mart = mart38)
out38 <- out38[!duplicated(out38$ensembl_gene_id),]
colnames(out38)[1] = "hg38"

hg19_38 <- full_join(out37, out38,  by = "ensembl_gene_id") %>% as.matrix()
hg19_38[is.na(hg19_38)] = ""


suppressMessages(MergedFPKMs$genes %<>% plyr::mapvalues(from = hg19_38[,"hg19"],
                                     to = hg19_38[,"ensembl_gene_id"]))
suppressMessages(MergedFPKMs$genes %<>% plyr::mapvalues(from = hg19_38[,"hg38"],
                                                        to = hg19_38[,"ensembl_gene_id"]))

suppressMessages(MergedFPKMs$genes %<>% plyr::mapvalues(from = hg19_38[,"hg19_transcript"],
                                       to = hg19_38[,"ensembl_gene_id"]))
fpkm_df <- fpkm_matrix %>% as.data.frame() %>% tibble::rownames_to_column("genes")
MergedFPKMsAll <- full_join(MergedFPKMs,fpkm_df,by = "genes")

suppressMessages(MergedFPKMsAll$genes %<>% plyr::mapvalues(from = hg19_38[,"ensembl_gene_id"],
                                         to = hg19_38[,"hg38"]))
MergedFPKMsAll %<>% filter(genes != "")

#=====remove duplicated genes with different ensembl_gene_id ====
table(duplicated(MergedFPKMsAll$genes))

genes <- MergedFPKMsAll$genes
dup <- genes[duplicated(genes)]
length(dup)
uniq <- genes[!(genes %in% dup)]
length(uniq)

MergedFPKMsAll_uniq <- MergedFPKMsAll[MergedFPKMsAll$genes %in% uniq,]
dim(MergedFPKMsAll_uniq)
MergedFPKMsAll_dup <- MergedFPKMsAll[MergedFPKMsAll$genes %in% dup,]
dim(MergedFPKMsAll_dup)
#----------- clean duplicated rows

MergedFPKMsAll_dup[is.na(MergedFPKMsAll_dup)] = 0

MergedFPKMsAll_dup_sum <- MergedFPKMsAll_dup %>% 
    split(f = MergedFPKMsAll_dup$genes) %>% 
    pbapply::pblapply(function(x) colSums(x[,-1])) %>% bind_rows()

dup_genes <- unique(MergedFPKMsAll_dup$genes)
table(dup_genes==dup)
MergedFPKMsAll_dup_sum <- cbind("genes" = dup_genes, MergedFPKMsAll_dup_sum)

MergedFPKMsAll_clean <- rbind(MergedFPKMsAll_uniq,MergedFPKMsAll_dup_sum)
MergedFPKMsAll_clean <- MergedFPKMsAll_clean[order(MergedFPKMsAll_clean$genes),]

test2 <- MergedFPKMsAll_clean[MergedFPKMsAll_clean$genes%in% markers,]
rownames(test2) = NULL
identical(test,test2[,1:(ncol(test2)-4)])

rownames(MergedFPKMsAll_clean) = NULL
MergedFPKMsAll_clean %<>% tibble::column_to_rownames("genes")
MergedFPKMsAll_clean[is.na(MergedFPKMsAll_clean)] = 0

write.csv(MergedFPKMsAll_clean,paste0(path,"MergedFPKMs.csv"),row.names = TRUE)

#=============

fpkm_matrix$genes <- plyr::mapvalues(rownames(fpkm_matrix),from = hg19_38[,"ensembl_gene_id"],
                                     to = hg19_38[,"hg38"])

Attr <- listAttributes(mart37)
Transcript <- Attr$name[grep("Transcript|transcript",Attr$description)]
List <- c(1:62)[-c(37,38,43,44,49,50,52)]

for(i in List) {
    print(i)
    test <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id',Transcript[i]),
                  filters = "ensembl_gene_id",values = "ENSG00000230999",
                  mart = mart37)

}

Transcript <- Transcript[-c(37,38,43,44,49,50,52)]

for(tran in Transcript){
    res <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id',tran),
          filters = "ensembl_gene_id",values = "ENSG00000230999",
          mart = mart37)
    print(paste(tran,"=",res[,3]))
}
