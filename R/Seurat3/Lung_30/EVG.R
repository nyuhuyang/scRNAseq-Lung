########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("dplyr","magrittr","data.table","pbapply","tibble","tidyr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
save.path <- "Yang/Lung_30/DE_analysis/F_EVGs/"
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

res_list <- vector(mode = "list",length = length(superfamily))
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
    res_list[[i]] = res
}

names(res_list) = superfamily
openxlsx::write.xlsx(res_list, file =  paste0(save.path,"Lung_30-EVGs-full.xlsx"),
                     colNames = TRUE,row.names = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

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

read.path = "Yang/Lung_30/DE_analysis/A_Sample_types/"
save.path = "Yang/Lung_30/DE_analysis/VolcanoPlots/"
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

DE_list <- c("Lung_30_A_57_celltypes=1_ALL CELLS_proximal_vs_distal+terminal",
             "Lung_30_A_113_celltypes=1_ALL CELLS_distal_vs_proximal+terminal",
             "Lung_30_A_225_celltypes=1_ALL CELLS_terminal_vs_proximal+distal",
             "Lung_30_A_281_celltypes=1_ALL CELLS_COPD_vs_distal")
for(i in seq_along(DE_list)) {
    DE_file = read.csv(paste0(read.path,DE_list[i],".csv"),row.names = 1)
    cell.type = sub("_vs_.*","",DE_list[i]) %>% gsub(".*_","",.)
    data = DE_file[DE_file$cluster %in% cell.type,]
    cut_off = "p_val_adj"
    cut_off_value = 0.05
    cut_off_logFC = 1
    top = 15
    sort.by = "p_val_adj"
    #cols = c("#2a71b2","#d2dae2","#ba2832")
    cols = c("dark green","light green","light grey","orange","red")
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
    data$change %<>% as.factor() %>% 
        factor(levels = c("Most Up-regulated",
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
    p<-ggplot(
        #设置数据
        data, 
        mapping = aes_string(x = "avg_logFC", 
                             y = paste0("log10_",cut_off[1]),
                             fill = "change"))+
        geom_point(alpha=alpha, size=size,color = "black", pch=21)
    # 辅助线
    p = p + geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8)
    p = p + geom_hline(yintercept = -log10(cut_off_value),lty=4,col="black",lwd=0.8) +
        
        # 坐标轴
        labs(x="log2(fold change)",
             y= paste("-log10 (",ifelse(cut_off[1] == "p_val_adj", "adjusted p-value","p-value"),")"))+
        
        # 图例
        theme(plot.title = element_text(hjust = 0.5), 
              axis.title=element_text(size=12),
              legend.position="right", 
              legend.title = element_blank(),
              legend.text = element_text(size = legend.size),
        )
    names(cols) = rev(c("Most Up-regulated",
                        "Up-regulated",
                        "Stable",
                        "Down-regulated",
                        "Most Down-regulated"))
    p = p +
        scale_fill_manual(values=cols[unique(data$change)])+ theme_bw()
    
    jpeg(paste0(save.path,"VolcanoPlots_",sub(".*celltypes=1_","",DE_list[i]),"_noLabel.jpeg"), family = "Arial",
         units="in", width=7, height=7,res=900)
    print(p)
    dev.off()
    p = p + ggrepel::geom_text_repel(data = data[c(Down_gene_index, Up_gene_index),], 
                                     aes(label = gene),
                                     size = 3,box.padding = unit(0.5, "lines"),
                                     point.padding = unit(0.8, "lines"), 
                                     segment.color = "black", 
                                     show.legend = FALSE)
    jpeg(paste0(save.path,"VolcanoPlots_",sub(".*celltypes=1_","",DE_list[i]),".jpeg"), family = "Arial",
         units="in", width=7, height=7,res=900)
    print(p)
    dev.off()
    svMisc::progress(i,length(DE_list))
}