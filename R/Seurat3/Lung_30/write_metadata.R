########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","magrittr","MAST","data.table","biomaRt","pbapply"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

save.path = "Yang/Lung_30/Cell_Phone_DB/"
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

#===== load data ===================
object = readRDS(file = "data/Lung_30_20200710.rds") 
print(unique(object@meta.data$conditions))
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
DefaultAssay(object) = "SCT"

#===== convert gene name to ENSEMBL===================

mart <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

ens <- rownames(object[["SCT"]]@data)
annotLookup <- getBM(mart=ensembl,
                     attributes=c("ensembl_gene_id","external_gene_name"),
                     filter="external_gene_name",
                     values=ens, uniqueRows=TRUE)
gene.keep = ens %in% annotLookup$external_gene_name
annotLookup = annotLookup[match(ens[gene.keep], annotLookup$external_gene_name),]


# For all
meta.data = cbind.data.frame(object@meta.data,
                             object@reductions$umap@cell.embeddings)
meta.data = meta.data[,c("UMAP_1","UMAP_2","SCT_snn_res.0.8",
                         "annotations3","orig.ident","conditions")]
colnames(meta.data) %<>% sub("orig.ident","samples",.)
colnames(meta.data) %<>% sub("conditions","regions",.)
colnames(meta.data) %<>% sub("annotations3","cell.types",.)
meta.data$SCT_snn_res.0.8 = as.numeric(as.character(meta.data$SCT_snn_res.0.8))
meta.data = meta.data[order(as.numeric(meta.data$SCT_snn_res.0.8)),]

write.csv(meta.data, paste0(save.path,"Cordinates_P_D_T_COPD.csv"),quote=F)

# For P, D, T separately =============
Idents(object) = "conditions"
#object$annotations3 %<>% gsub("-.*","",.) %>% gsub("[0-9]+","",.)
conditions <- c("proximal","distal","terminal","COPD")
for (con in  conditions){
        message(paste("write meta.data for", con))
        sub_object <- subset(object,idents = con)
        colnames(sub_object@meta.data) %<>% sub("annotations3","cell.types",.)
        
        #meta
        meta.data = cbind.data.frame(sub_object@meta.data,
                                     sub_object@reductions$umap@cell.embeddings)
        meta.data = meta.data[,c("UMAP_1","UMAP_2","SCT_snn_res.0.8",
                                 "cell.types","orig.ident","conditions")]
        colnames(meta.data) %<>% sub("orig.ident","samples",.)
        colnames(meta.data) %<>% sub("conditions","regions",.)

        meta.data$SCT_snn_res.0.8 = as.numeric(as.character(meta.data$SCT_snn_res.0.8))
        
        meta.data = meta.data[order(as.numeric(meta.data$SCT_snn_res.0.8)),]
        #write.table(meta.data, paste0(save.path,"Lung_30-",con,"_meta.data.txt"), sep='\t', quote=F)
        
        meta_data <- cbind("Cell" = rownames(sub_object@meta.data), sub_object@meta.data[,"cell.types", drop=F])   #####  cluster is the user’s specific cluster column
        write.table(meta_data, paste0(save.path,"Lung_30-",con,"_meta.data_cellphone.txt"), sep='\t', quote=F, row.names=F)
        
        message(paste("write counts for", con))
        counts = as.data.table(sub_object@assays$SCT@data,keep.rownames=TRUE)
        rownames(counts) = rownames(sub_object@assays$SCT@data)
        colnames(counts)[1] = "Gene"
        counts = counts[gene.keep,]
        print(table(counts$Gene == annotLookup$external_gene_name))
        counts$Gene = annotLookup$ensembl_gene_id
        
        fwrite(counts,
               file = paste0(save.path,"Lung_30-",con,"_counts.txt"),
               sep = "\t",row.names = FALSE, col.names = TRUE)
        Progress(which(conditions %in% con) , length(conditions))
        }

object_list <- SplitObject(object, split.by = "conditions")
counts_list <- lapply(object_list, function(x) x@assays$SCT@counts)
dt_list <- lapply(counts_list, function(x) as.data.table(x,keep.rownames=TRUE))
(conditions <- names(counts_list))

for(i in seq_along(conditions)){
        fwrite(dt_list[[i]],file = paste0(save.path,"Lung_30-",conditions[i],"_counts.txt"),
               sep = "\t",row.names = FALSE, col.names = TRUE)
        Progress(i, length(conditions))
}


#============= conditions===========
# Columns for each cell type separately (as many columns as cell types) 
# to show average expression in each cell type,
conditions <- c("proximal","distal","terminal","COPD")
Idents(object) = "conditions"
Idents(object) %<>% factor(levels = conditions)
conditions_exp <- AverageExpression(object = object, 
                          features = VariableFeatures(object),
                          assays = "SCT")
conditions.pct <- pbsapply(conditions, function(x)  {
                                thresh.min <- 0
                                features = VariableFeatures(object)
                                cells.1 = colnames(object)[object$conditions %in% x]
                                pct.1 <- round(
                                        x = rowSums(x = object@assays$SCT[features, cells.1] > thresh.min) /
                                                length(x = cells.1),
                                        digits = 3)
                                }
                            )
colnames(conditions.pct) %<>% paste0(".pct")
conditions_exp$SCT %<>% cbind(conditions.pct)
conditions_exp = conditions_exp$SCT[,c("proximal","proximal.pct","distal","distal.pct","terminal","terminal.pct","COPD","COPD.pct")]
#  + 4 extra columns for average expression in P (all cells - 1 column),
# D (all cells - 1 column), T (all cells - 1 column), COPD (all cells - 1 column). 
df_abbr <- readxl::read_excel("doc/Annotations/Cell type abbreviation.xlsx",
                              sheet = "abbreviation_20201116")

table(df_abbr$Abbreviation %in% object$annotations3)
object$cell_types <- plyr::mapvalues(object$annotations3,
                                            from = df_abbr$Abbreviation,
                                            to = df_abbr$`Revised abbreviations`)

Idents(object) = "cell_types"
object %<>% sortIdent()
cell_types = sort(unique(object$cell_types))

cell_types_exp <- AverageExpression(object = object, 
                                 features = VariableFeatures(object),
                                 assays = "SCT")
cell_types.pct <- pbsapply(cell_types, function(x)  {
                                thresh.min <- 0
                                features = VariableFeatures(object)
                                cells.1 = colnames(object)[object$cell_types %in% x]
                                pct.1 <- round(
                                        x = rowSums(x = object@assays$SCT[features, cells.1] > thresh.min) /
                                                length(x = cells.1),
                                        digits = 3)
}
)
colnames(cell_types.pct) %<>% paste0(".pct")
cell_types_exp$SCT %<>% cbind(cell_types.pct)
cell_types_exp = cell_types_exp$SCT[,sort(colnames(cell_types_exp$SCT))]
Expression <- cbind(cell_types_exp, conditions_exp)

# number of genes expressed in each cell type vs conditions.
Orig.ident = c("UNC-48-P","UNC-55-P","UNC-66-P","UNC-69-P","UNC-71-P","UNC-48-D","UNC-55-D",
               "UNC-66-D","UNC-69-D","UNC-71-D","UNC-44-D","UNC-54-D","UNC-57-D","UNC-67-D",
               "UNC-70-D","CU-12-D","CU-12-D-R","VU-27-D","UNC-48-T","UNC-55-T","UNC-66-T",
               "UNC-69-T","UNC-71-T","UNC-44-T","CU-12-T","UNC-51-D","UNC-52-D","UNC-61-D",
               "VU19-D","VU-37-D")
nFeature_SCT_all <- pbsapply(cell_types, function(x){
        round(mean(object@meta.data[object$cell_types %in% x, "nFeature_SCT"]),digits = 0)
        })
nFeature_SCT_samples <- pbsapply(Orig.ident, function(y){
        pbsapply(cell_types, function(x){
                round(mean(object@meta.data[(object$cell_types %in% x) & (object$orig.ident %in% y),
                              "nFeature_SCT"]),digits = 0)
                })
        })
nFeature_SCT_samples[is.nan(nFeature_SCT_samples)] = 0
avg.number_of_genes_expressed = cbind("All"=nFeature_SCT_all,
                                      nFeature_SCT_samples[,Orig.ident])
n_cells = data.frame("All" = unclass(table(object$cell_types)))
orig.ident_cells = as.data.frame(unclass(table(object$cell_types, object$orig.ident)))

# number of cells in each cell type vs orig.ident
total.number_of_cells = cbind(n_cells,orig.ident_cells[,Orig.ident] )
output = list("avg.number_of_genes_expressed"= avg.number_of_genes_expressed,
              "total.number_of_cells"=total.number_of_cells)
#=====
data.table::fwrite(Expression,file = paste0(path,"Lung_30-Cell.types_expression.txt"),
                   sep = "\t",row.names = TRUE, col.names = TRUE)
openxlsx::write.xlsx(output, 
                     file =  paste0(path,"Lung_30-avg.number_of_genes_expressed.xlsx"),
                     colNames = TRUE,row.names = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))
