########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#conda activate r4.0.3 
# devtools::install_github('satijalab/seurat-data') #3.1.5.9900
# remotes::install_github("mojaveazure/seurat-disk")
invisible(lapply(c("Seurat","dplyr","SeuratData","SeuratDisk","cowplot","magrittr","tidyr","tibble"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
save_path <- "Yang/Lung_30/Deconvolution/"
if(!dir.exists(save_path))dir.create(save_path, recursive = T)

# 5.1.1 load data
# Rename ident
object = readRDS(file = "data/Lung_SCT_30_20200710.rds")
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
Idents(object) = "cell_types"
object %<>% subset(idents = "Un", invert = T)
df_samples <- readxl::read_excel("doc/Annotations/Cell type abbreviation - refs for deconvolution.xlsx")
table(df_samples$cell_types %in% object@meta.data[,"cell_types"])

groups = c("cell_types.select","cell_types","major.cell_types","RS.group","cell.family","super.family")
for(g in groups){
        object[[g]] = plyr::mapvalues(object@meta.data[,"cell_types"],
                                      from = df_samples$cell_types,
                                      to = df_samples[[g]])
}

for(g in groups){
        freq_table <- prop.table(x = table(object@meta.data[[g]],object$orig.ident),
                                 margin = 2) %>% as.data.frame.matrix
        write.table(t(freq_table), paste0(save_path,"psudobulk_",g,"_percent.txt"),sep = "\t",quote = FALSE)
}
# remove "---" column
# generate mixture
Idents(object) = "orig.ident"
expression = AggregateExpression(object,assay = "SCT")
expression = expression$SCT
tpm.mat = t(t(expression)*1e4 /colSums(expression) )
write.table(tpm.mat,file = paste0(save_path,"psudobulk_tpm_Lung30_orig.ident.txt"),sep = "\t",quote = FALSE)
# add "GeneSymbol" to the first line.

# convert to h5ad
object[["pca"]] = NULL
SaveH5Seurat(object, filename = "data/lung_30.h5Seurat")
Convert("data/lung_30.h5Seurat", dest = "h5ad")


meta_color = table(object$cell_types, object$cell_types.colors) %>% t %>% as.data.frame.matrix()
((meta_color > 0)*1) %>% rowSums()
UMAPPlot(object,reduction = NULL, cols  =ExtractMetaColor(object)) + NoLegend()

meta.data = object@meta.data
meta.data = meta.data[!duplicated(meta.data$cell_types.colors),]

freq_table = freq_table *100
Cell_proportion <- readxl::read_excel(paste0("doc/","Cell proportion - single cell data.xlsx"))
table(Cell_proportion$cell.types %in% rownames(freq_table))
table(colnames(Cell_proportion) %in% colnames(freq_table))

freq_table = freq_table[Cell_proportion$cell.types,colnames(Cell_proportion)[-1]]
write.table(freq_table, paste0(path,"psudobulk_cell_type_percentage.txt"),sep = "\t",quote = FALSE)
freq_table$cell_types = rownames(freq_table)

freq_table %<>% gather("sample","percentage",-cell_types)
freq_table$cell_types.color = plyr::mapvalues(freq_table$cell_types, from  = meta.data$cell_types,
                                              to = meta.data$cell_types.colors)
ggOut = ggplot(freq_table, aes_string(x = "sample",
                                      y =    "percentage",
                                      fill = "cell_types")) +
        geom_col() + 
        ggtitle("ground truth")+
        xlab("sample")+
        ylab("Cell Proportion (%)")+#NoLegend()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size =10),
              plot.title = element_text(hjust = 0.5))+
        scale_fill_manual(values = ExtractMetaColor(object,group.by = "cell_types"))
ggOut

cross.entropy <- function(p, phat){
        x <- 0
        for (i in 1:length(p)){
                if(phat[i] == 0) next
                x <- x + (p[i] * log(phat[i]))
        }
        return(-x)
}
cross.entropy(freq_table[freq_table$sample %in% "CU-12-D","percentage"],
              freq_table[freq_table$sample %in% "CU-12-D","percentage"]
)


# generate reference
for(g in groups){
        print(g)
        Idents(object) = g
        if("---" %in% unique(Idents(object))){
                sub_object = subset(object, idents = "---", invert = T)
        } else sub_object = object
        
        expression = AggregateExpression(sub_object,assay = "SCT",group.by = g)
        expression = expression$SCT
        tpm.mat = t(t(expression)*1e6 /colSums(expression) )
        write.table(tpm.mat,file = paste0(save_path,"psudobulk_tpm_Lung30_",g,"_reference.txt"),sep = "\t",quote = FALSE)
        }


# add "gene" to the first line

# generate single-cell reference
counts = object[["SCT"]]@counts
cell_types <- sort(unique(object$cell_types))
meta.data = object@meta.data
select_barcode  <- pbapply::pbsapply(cell_types, function(cell_type){
        sub_meta.data = meta.data %>% filter(cell_types %in% cell_type)
        sample(rownames(sub_meta.data), size = min(nrow(sub_meta.data),100))
        
})
select_barcode %<>% unlist()
sub_object = subset(object,cells = select_barcode)

counts = as.matrix(sub_object[["SCT"]]@counts)
sub_meta.data =sub_object@meta.data
colnames(counts) = sub_object$cell_types
write.table(counts,file = paste0(path,"single.cell_Lung30_cell.type_reference.txt"),sep = "\t",quote = FALSE)

# add "GeneSymbol" to the first line.

# CIBERSORT
CIBERSORT <- read.csv("./output/CIBERSORT.csv",row.names = 1)
CIBERSORT <- CIBERSORT[,-((ncol(CIBERSORT)-2):ncol(CIBERSORT))]
t(CIBERSORT)[1:5,]
# 5.1.2 EpiDISH
library(EpiDISH)
MCL.expression <- read.delim2("./output/MCL_expression.txt",row.names = 1)
refer <- read.delim2("./data/refer.txt",row.names = 1)
MCL.expression <- data.matrix(MCL.expression) # convert dataframe to matrix
refer <- data.matrix(refer) # convert dataframe to matrix
head(sapply(refer,class))
head(sapply(MCL.expression,class))

out.l <- epidish(MCL.expression, refer, method = "RPC", maxit = 100) 
results_RPC <- out.l$estF
t(results_RPC)[1:5,]

out.l <- epidish(MCL.expression, refer, method = "CBS") 
results_CBS <- out.l$estF
t(results_CBS)[1:5,]

#out.l <- epidish(avdata.m=MCL.expression, ref.m=refer, method = "CP",constraint = "equality") 
#results_CP <- out.l$estF

# 5.1.3 DeconRNASeq
library(DeconRNASeq)
MCL.expression <- read.delim2("./output/MCL_expression.txt",row.names = 1)
refer <- read.delim2("./data/refer.txt",row.names = 1)
MCL.expression <- as.data.frame(data.matrix(MCL.expression)) # convert to dataframe 
refer <-  as.data.frame(data.matrix(refer)) # convert to dataframe 
head(sapply(refer,class)) 
head(sapply(MCL.expression,class))

results_Decon <- DeconRNASeq(MCL.expression, refer, NULL, checksig=FALSE, 
            known.prop = FALSE, use.scale = FALSE, fig = FALSE)
rownames(results_Decon$out.all) <- colnames(MCL.expression)
results_Decon <- results_Decon$out.all
t(results_Decon)[1:5,]

# 5.1.4 xCell
library(xCell)
par(mar=c(14,3,3,3))
expression <- read.table("output/20210726/psudobulk_tpm_Lung30_orig.ident_mixture.txt",header=TRUE,row.names=1, as.is=TRUE)
refer <- read.table("output/20210726/psudobulk_tpm_Lung30_cell.type_reference.txt",header=TRUE,row.names=1, as.is=TRUE)
refer <- data.matrix(refer) # convert dataframe to matrix
head(sapply(refer,class))

# convert dataframe column to list r
genes <- rownames(refer)
refer.list <- lapply(1:ncol(refer), function(x) refer[,x])
refer.list <- setNames(refer.list, colnames(refer))
lapply(refer.list,head)
refer.list <- lapply(refer.list, function(x) x[order(x,decreasing = T)])
lapply(refer.list,head)
refer.list <- lapply(refer.list, function(x) names(x[1:300]))
lapply(refer.list,length)
genes <-unique(unlist(refer.list))
length(genes)
expression = expression[complete.cases(expression),]
results_xCell <- xCellAnalysis.1(expr=expression,signatures = refer.list)
# divide each row of a matrix by elements of a vector in R
results_xCell <- t(t(results_xCell)/colSums(results_xCell))
colSums(results_xCell)
results_xCell1 <- as.data.frame(results_xCell)
results_xCell1$cell_types = gsub("\\.","-",rownames(results_xCell1))
results_xCell1 %<>% gather("sample","percentage",-cell_types)

cell_types_color = ExtractMetaColor(object,group.by = "cell_types")

results_xCell1$cell_types.color = plyr::mapvalues(results_xCell1$cell_types, 
                                                  from  = names(cell_types_color),
                                              to = cell_types_color)
ggOut = ggplot(results_xCell1, aes_string(x = "sample",
                                      y =    "percentage",
                                      fill = "cell_types")) +
        geom_col() + 
        ggtitle("xcell deconvolution")+
        xlab("sample")+
        ylab("Cell Proportion (%)")+#NoLegend()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              plot.title = element_text(hjust = 0.5))+
        scale_fill_manual(values = cell_types_color)
ggOut


ggOut

# 5.1.5 DeMixT
library("DeMixT")
MCL.expression <- read.delim2("./output/MCL_expression.txt",row.names = 1)
refer <- read.delim2("./data/refer.txt",row.names = 1)
MCL.expression <- data.matrix(MCL.expression) # convert dataframe to matrix
refer <- data.matrix(refer) # convert dataframe to matrix
head(sapply(refer,class))
head(sapply(MCL.expression,class))

#Error in cbind(data.comp1, data.Y) : number of rows of matrices must match (see arg 2)
CommonGenes <- intersect(rownames(MCL.expression),rownames(refer))
MCL.expression <- MCL.expression[CommonGenes,]
refer <- refer[CommonGenes,]
dim(MCL.expression)
dim(refer)
gene.use <- apply(refer,1,sd)
gene.use <- gene.use[order(gene.use,decreasing = T)]
plot(gene.use)
refer.highsd <- refer[gene.use[1:100],]
MCL.expression.highsd <- MCL.expression[gene.use[1:100],]
dim(MCL.expression.highsd)
dim(refer.highsd)

results_DeMixT <- DeMixT(data.Y = MCL.expression.highsd, data.comp1 = refer.highsd,
                         if.filter = FALSE) 
results_DeMixT <- DeMixT.S1(data.Y = MCL.expression, data.comp1 = refer, if.filter = FALSE) 

# 5.2 Summarize all deconvolution results
rownames(CIBERSORT) <- paste0(rownames(CIBERSORT),".CIBERSORT")
rownames(results_RPC) <- paste0(rownames(results_RPC),".epidish_RPC")
rownames(results_CBS) <- paste0(rownames(results_CBS),".epidish_CBS")
rownames(results_Decon) <- paste0(rownames(results_Decon),".DeconRNASeq")

results <- data.frame(t(CIBERSORT),t(results_RPC),t(results_CBS),t(results_Decon))
results_all <- full_join(rownames_to_column(results), 
                         rownames_to_column(Seurat_freq), by = "rowname")
results_all$Cell_Types <- gsub('\\..*', '',results_all$rowname)
results_all <- results_all[,-1]
results_all[is.na(results_all)] <- 0
results_short <- aggregate(. ~ Cell_Types, data = results_all, FUN=sum)
rownames(results_short) <- results_short$Cell_Types
head(results_short)
#results_short <- results_short[,-1]
#boxplot(t(results_short))

results_short2 <- melt(results_short, value.name = "percentage")
results_short2$variable <- as.character(results_short2$variable)

Catalog <- matrix(unlist(strsplit(results_short2$variable,".",fixed = T)),
               ncol=2,byrow = T)
colnames(Catalog) <- c("patient_vs_normal","software")
results_short2 <- cbind.data.frame(results_short2,Catalog)
head(results_short2)
results_short2$Cell_Types <- as.factor(results_short2$Cell_Types)
ggplot(results_short2,aes(x = Cell_Types, y = percentage,
                                  group=software))+
        geom_col(aes(fill = Cell_Types))+
        facet_grid(software~patient_vs_normal)+
        scale_y_continuous(labels = scales::percent_format())+
        geom_text(aes(label = scales::percent(percentage)), vjust = -.5)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
