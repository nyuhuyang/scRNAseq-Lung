########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("dplyr","magrittr","data.table","pbapply","tibble","tidyr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
read.path = "Yang/Lung_30/DE_analysis/A_Sample_types/"
#read.path = "Yang/Lung_30/DE_analysis/D_age/"

save.path = "Yang/Lung_30/DE_analysis/F_EVGs_allCells/"
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

#=======
Expression = data.table::fread(file = "output/20201215/Lung_30-Cell.types_expression.txt",
                   sep = "\t") %>% as.data.frame()
Expression %<>% column_to_rownames(var = "V1")
Expression[1:4,1:4]

anno <- readxl::read_excel("doc/Annotations/Cell type abbreviation.xlsx")

df <- readxl::read_excel("doc/Chord diagram cell order - updated abbreviations 12-14-20.xlsx",col_names = T)
superfamily <- c("Epithelial","Structural","Immune")
EVG <- vector(mode = "list",length = length(superfamily))
for(i in seq_along(superfamily)){
    selected.cell.types <- df$cell.types[df$Cell.group %in% superfamily[i]] %>% sort
    gene_list = pbsapply(selected.cell.types,function(ident.1){
        ident.2 = selected.cell.types[!(selected.cell.types %in% ident.1)]
        keep1 = Expression[,paste0(ident.1,".pct")] >= 1/3
        keep2 = Expression[,ident.1]/rowMeans(Expression[,ident.2]) >= 2
       rownames(Expression)[keep1 & keep2]
    })
    gene_list %<>% list2df
    gene_list[is.na(gene_list)] = ""
    EVG[[i]] = gene_list
}
names(EVG) = superfamily
openxlsx::write.xlsx(EVG, file =  paste0(save.path,"Lung_30-EVGs.xlsx"),
                     colNames = TRUE,row.names = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))


#===== generate FC===================
EVG_df <- pblapply(EVG,function(x) {
    temp = gather(as.data.frame(x),cell.type,genes)
    temp[temp$genes != "",]
    }) %>% bind_rows

EVGs_list <- vector(mode = "list",length = length(superfamily))
cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix("", n-nrow(x), ncol(x))))) 
}

for(i in seq_along(superfamily)){
    selected.cell.types <- df$cell.types[df$Cell.group %in% superfamily[i]] %>% sort
    genes.de = pblapply(selected.cell.types,function(ident.1){
        ident.2 = selected.cell.types[!(selected.cell.types %in% ident.1)]
        genes  = EVG_df %>% filter(cell.type == ident.1) %>% dplyr::select(genes) %>% pull
        
        keep1 = Expression[,paste0(ident.1,".pct")] >= 1/3
        keep2 = Expression[,ident.1]/rowMeans(Expression[,ident.2]) >= 2
        rownames(Expression)[keep1 & keep2]
        de = data.frame(gene  = genes,
                        exp.1 = Expression[genes,ident.1],
                        exp.2 = rowMeans(Expression[genes,ident.2]),
                        pct.1 = Expression[genes,paste0(ident.1,".pct")],
                        pct.2 = rowMeans(Expression[genes,paste0(ident.2,".pct")]),
                        FC = Expression[genes,ident.1]/rowMeans(Expression[genes,ident.2])
        )
        de %<>% arrange(desc(FC))
        colnames(de) %<>% paste0(ident.1,"_",.)
        return(de)
    })
    res <- do.call(cbind.fill, genes.de)
    rownames(res) = NULL
    EVGs_list[[i]] = res
}

names(EVGs_list) = superfamily
openxlsx::write.xlsx(EVGs_list, file =  paste0(save.path,"Lung_30-EVGs-full.xlsx"),
                     colNames = TRUE,row.names = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))
#========= expression data for 577 GTEx samples ==============
library(data.table)
EVGs_df <- do.call(cbind.fill, EVGs_list)
EVGs_sorted = EVGs_df[,grep("gene",colnames(EVGs_df))]
colnames(EVGs_sorted) %<>% gsub("_gene$","",.)
cell.types = colnames(EVGs_sorted)
df_TPM = fread("data/RNA-seq/GTEx-Lung-tpm.csv",header = T)
if(df_TPM$V1[1] == "Age.Bracket") df_TPM = df_TPM[-1,]
TPM <- column_to_rownames(df_TPM, var = "V1") %>% sapply(as.numeric) %>%
    as.data.frame
TPM = cbind("V1"=df_TPM$V1,TPM)
rownames(TPM) = TPM$V1

temp = readxl::read_excel("Yang/GTEx/Cell type EVG genes - 577 GTEx samples ordered.xlsx",sheet = "AT1",)
colnames(temp)[1] = "V1"
exp_list <- list()
for(i in seq_along(cell.types)) {
    genes <- EVGs_sorted[,cell.types[i]]
    genes = genes[genes != ""]
    genes = genes[!is.na(genes)]
    exp <- TPM[genes,colnames(temp)]
    genes = genes[genes %in% exp$V1]
    exp = exp[match(genes, exp$V1),]
    exp_list[[i]] = exp#rbind.data.frame(temp,exp)
    colnames(exp_list[[i]])[1] = ""
    Progress(i, length(cell.types))
}

names(exp_list) = cell.types
openxlsx::write.xlsx(exp_list, file =  "Yang/GTEx/Cell type EVG genes - 577 GTEx samples ordered~.xlsx",
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

# =============Volcano plots =============
#Include all genes

#Formatting:

#“Most” Up-regulated (log2 FC >1) in the first group vs compared group: red dots, and include their gene names.
#“Most” down-regulated (log2 FC < -1): green dots, and include their names.
#Remaining up-regulated (p adj <0.05) - orange  
#Remaining down-regulated (p adj <0.05) - light green 
#Remaining (p adj <0.05) - light grey
#Borders of dots - black, thickness 0.2

#Aspect ratio: the figures will have approximate size in the figure: height 1.2 width 1.2.

read.path = "Yang/Lung_30/DE_analysis/F_EVGs_allCells/"
save.path = "Yang/Lung_30/DE_analysis/VolcanoPlots/F_EVGs/"
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

DE_list <- list.files(path = read.path,pattern = "\\.csv")

# read EVGs
superfamily <- c("Epithelial","Structural","Immune")
EVGs_df <- lapply(superfamily, function(s) {
    tmp = readxl::read_excel(paste0(read.path,"Lung_30-EVGs-full.xlsx"), sheet = s)
    tmp = tmp[,-1]
    tmp = tmp[,grep("_gene",colnames(tmp),value = T)]
    colnames(tmp) %<>% gsub("_gene","",.)
    tmp
}) %>% do.call(cbind.fill,.)


for(i in seq_along(DE_list)) {
    DE_file <- tryCatch(expr = read.csv(paste0(read.path,DE_list[i]),row.names = 1),
                error = function(cond) {
        return(cond$message)
    })
    
    (cell.type <- sub("_vs_.*","",DE_list[i]) %>% 
            stringi::stri_split_fixed(pattern = "_",simplify = T) %>% .[5])
    ident1 = sub("_vs_.*","",DE_list[i]) %>% sub(".*_","",.)
    ident2 = sub(".*_vs_","",DE_list[i]) %>% sub("_.*","",.) %>% sub(".csv","",.)
    
    data = DE_file[DE_file$cluster %in% ident1,]
    
    EVGs = EVGs_df[,cell.type]
    EVGs = EVGs[!(EVGs == "" | is.na(EVGs))]
    
    cut_off = "p_val_adj"
    cut_off_value = 0.05
    cut_off_logFC = 1
    top = 15
    sort.by = "p_val_adj"
    cols = rev(c("dark green","light green","light grey","orange","red","black"))
    names(cols) = c("Enriched Variable Genes(EVGs)",
                    "Most Up-regulated",
                    "Up-regulated",
                    "Stable",
                    "Down-regulated",
                    "Most Down-regulated")
    alpha=0.9
    size=2
    
    data[,paste0("log10_",cut_off[1])] = -log10(data[,cut_off[1]])
    data[,"change"] = "Stable"
    data[data$p_val_adj <= cut_off_value &
             data$avg_logFC > 0,"change"] = "Up-regulated"
    data[data$p_val_adj <= cut_off_value &
             data$avg_logFC < 0,"change"] = "Down-regulated"
    
    data[data$avg_logFC > cut_off_logFC &
             data$p_val_adj < cut_off_value,"change"] = "Most Up-regulated"
    data[data$avg_logFC < -cut_off_logFC &
             data$p_val_adj < cut_off_value,"change"] = "Most Down-regulated"
    
    # now show EVGs if Stable
    #EVGs = EVGs[!(EVGs %in% data[data$change %in% "Stable","gene"])]
    if(length(EVGs)>0) {
        EVGs = EVGs[1:min(length(EVGs),top)]
        data[EVGs,"change"] = "Enriched Variable Genes(EVGs)"
    }
    data$change %<>% as.factor() %>% 
        factor(levels = c("Enriched Variable Genes(EVGs)",
                          "Most Up-regulated",
                          "Up-regulated",
                          "Stable",
                          "Down-regulated",
                          "Most Down-regulated"))
    colnames(data)[grep("cluster",colnames(data))]="cluster"
    # 将需要标记的基因放置在单独的数组
    Up <- data[data$change %in% "Most Up-regulated",]
    Down <- data[data$change %in% "Most Down-regulated",]
    if(sort.by == "p_val_adj") {
        Up_gene_index <- rownames(Up)[Up[,sort.by] <= tail(head(sort(Up[,sort.by],decreasing = F),top),1)]
        Down_gene_index <- rownames(Down)[Down[,sort.by] <= tail(head(sort(Down[,sort.by],decreasing = F),top),1)]
    }
    if(sort.by == "avg_logFC") {
        Up_gene_index <- rownames(Up)[Up[,sort.by] >= tail(head(sort(Up[,sort.by],decreasing = T),top),1)]
        Down_gene_index <- rownames(Down)[Down[,sort.by] <= tail(head(sort(Down[,sort.by],decreasing = F),top),1)]
    }
    # If too many gene with adj_p = 0
    if(length(Up_gene_index) > top) Up_gene_index %<>% sample(top)
    if(length(Down_gene_index) > top) Down_gene_index %<>% sample(top)
            
    p <- ggplot(
        #设置数据
        data, 
        mapping = aes_string(x = "avg_logFC", 
                             y = paste0("log10_",cut_off[1]),
                             fill = "change"))+
        geom_point(alpha=alpha, size=size,color = "black", pch=21)+
        ggtitle(label = paste(ident1, "vs.", ident2))
    # 辅助线
    p = p + geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8)
    p = p + geom_hline(yintercept = -log10(cut_off_value),lty=4,col="black",lwd=0.8) +
        
        # 坐标轴
        labs(x="log2(fold change)",
             y= paste("-log10 (",ifelse(cut_off[1] == "p_val_adj", "adjusted p-value","p-value"),")"))+ 
        theme_bw()+
        
        # 图例
        theme(plot.title = element_text(hjust = 0.5), 
              axis.title=element_text(size=12),
              legend.position="right", 
              legend.title = element_blank(),
              legend.text = element_text(size = 10),
        )
    p = p +
        scale_fill_manual(values=cols[unique(data$change)])
    
    p = p + ggrepel::geom_text_repel(data = data[c(Down_gene_index, Up_gene_index,EVGs),], 
                                     aes(label = gene),
                                     size = 3,box.padding = unit(0.5, "lines"),
                                     point.padding = unit(0.8, "lines"), 
                                     segment.color = "black", 
                                     show.legend = FALSE)
    fileName = stringi::stri_split_fixed(basename(DE_list[i]),pattern = "_",simplify = T,n = 5)[5] %>% sub("\\.csv","",.)
    jpeg(paste0(save.path,"VolcanoPlots_",fileName,".jpeg"), family = "Arial",
         units="in", width=7, height=7,res=900)
    print(p)
    dev.off()
    svMisc::progress(i,length(DE_list))
}

#============
(DE_list <- list.files(path = read.path,pattern = "\\.csv") %>% grep("ALL CELLS",.,value = T))
(DE_list <- list.files(path = read.path,pattern = "\\.csv") %>% grep("ALL ",.,value = T))

EVGs_long <- gather(as.data.frame(EVGs_df), "cell.type","gene") %>%
    filter(gene != "")
EVGs_full <- data.frame("gene" = sort(unique(EVGs_long$gene)))
cell.types <- unique(EVGs_long$cell.type)
for(i in seq_along(cell.types)) {
    tmp = EVGs_long[EVGs_long$cell.type %in% cell.types[i],]
    EVGs_full %<>% left_join(tmp,by = "gene")
}
EVGs_full %<>% as.matrix()
EVGs <- apply(EVGs_full,1, function(x) as.character(na.exclude(x))) %>%
    do.call(cbind.fill,.) %>% t %>% as.data.frame()
colnames(EVGs) = 1:(ncol(EVGs))-1
colnames(EVGs) %<>% paste0("EVGs_",.)
colnames(EVGs)[1] = "gene"
for(i in seq_along(DE_list)) {
    DE_file <- tryCatch(expr = read.csv(paste0(read.path,DE_list[i]),row.names = 1),
                        error = function(cond) {
                            return(cond$message)
                        })
    DE_file %<>% full_join(EVGs,by = "gene") %>% as.matrix()
    
    DE_file[is.na(DE_file)]=""
    write.csv(DE_file,paste0(save.path,DE_list[i]))
    svMisc::progress(i,length(DE_list))
    }
