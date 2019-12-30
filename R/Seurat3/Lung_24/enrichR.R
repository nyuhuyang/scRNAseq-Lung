########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(ggplot2)
library(enrichR)
library(tibble)
library(ggpubr)
library(ggsci)
source("../R/Seurat3_functions.R")
path <- paste0("Yang/proximal_distal_terminal/Non-Integration/GSEA")
if(!dir.exists(path)) dir.create(path, recursive = T)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
DEG_path <- "Yang/proximal_distal_terminal/Non-Integration/DEGs/cell_types/"

# read csv files ===============
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
df_cell_types$`Cell types` %<>%  gsub("/","-",.) #NK-T

CSV_list = list.files(path = DEG_path, pattern = "Lung_24-FC0_markers_")
cell_types <- gsub(".*_","",CSV_list) %>% gsub(".csv","",.) %>% sort
CSV_files <- list()
for (i in seq_along(cell_types)) {
        CSV_list[i] = list.files(path = DEG_path, pattern = cell_types[i])
        CSV_list
        CSV_files[[i]] = read.csv(paste0(DEG_path,CSV_list[i]), row.names = 1,stringsAsFactors = F)
        CSV_files[[i]]$cluster = gsub(".*_","",CSV_list[i]) %>% gsub(".csv","",.)
        CSV_files[[i]]$gene = rownames(CSV_files[[i]])
}
CSV_files <- bind_rows(CSV_files)

table(df_cell_types$`Cell types` %in% CSV_files$cluster)

write.csv(CSV_files,file = paste0(DEG_path,"Lung_24-FC0_cell_types.csv"))
CSV_files= read.csv(file= paste0(DEG_path,"Lung_24-FC0_cell_types.csv"))
CSV_files$cluster %<>% plyr::mapvalues(from = df_cell_types$`Cell types`,
                                       to = df_cell_types$Abbreviation)

# load pathway database ===============
dbs_table <- readxl::read_excel("Yang/proximal_distal_terminal/Non-Integration/GSEA/cell types/Enrichr libraries.xlsx",
                                sheet = "Sheet2")
dbs <- dbs_table$`Gene-set Library`

CSV_files = CSV_files[CSV_files$p_val_adj <0.05,]
CSV_files = CSV_files[CSV_files$avg_logFC >= 1,]

#' FgseaDotPlot generate Dot plot using findmarker results based on FGSEA
#' @param stats findmarker results
#' @param pathways pathway database list
#' @param Rowv determines if and how the row dendrogram should be reordered. 
#' By default, NULL or FALSE, then no dendrogram is computed and no reordering is done.
#' If a vector of integers, then dendrogram is computed and reordered based on the order of the vector.
#' @param Colv 	determines if and how the column dendrogram should be reordered.
#' Has the options as the Rowv argument above.
#' @param title add to title names
#' @param padj padj cut off
#' @param pval pval cut off
#' @param order.yaxis.by c(1,"pval") means order y axis by pval in cluster 1
#' @param order specify order of x axis
#' @param nperm fgsea param, default 1000
#' @param do.return return fgsea data frame
#' @param return.raw return fgsea raw data
#' @param ... ggballoonplot param
#' @example FgseaDotPlot(stats=res, pathways=hallmark,title = "each B_MCL clusters")
EnrichDotPlot <- function(stats=CSV_files, pathways=NULL,
                         size = "-log10(pval)", Rowv = NULL,Colv = NULL,
                         font.ytickslab = 15,
                         fill = "NES",title="",order.row.by = c(1,"Adjusted.P.value"),
                         order.col = NULL,decreasing = T,
                         pathway.name = "Hallmark",padj = 0.25, pval=0.05,
                         do.return = F,return.raw = F,font.main = 18,
                         verbose=T,save_path = NULL, width=10, height=7,hjust=0.5,...){
        
        clusters =  stats$cluster %>% as.character %>% unique %>% sort
        message("Calculate fgsea for each cluster.")
        enrichedRes <- list()
        for(i in seq_along(clusters)){
                geneRank = stats[stats$cluster == clusters[i],]
                geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")]  %>% tibble::deframe()
                enrichedRes[[i]] <- enrichr(names(geneRank), dbs)
                for(k in seq_along(enrichedRes[[i]])) {
                        if(nrow(enrichedRes[[i]][[k]]) > 0 ) enrichedRes[[i]][[k]][,"library"] = names(enrichedRes[[i]][k])
                }
                enrichedRes[[i]] = bind_rows(enrichedRes[[i]])
                if(clusters[i] == order.row.by[1]) {
                        order.row = enrichedRes[[i]][order(enrichedRes[[i]][,order.row.by[2]],
                                                        decreasing = decreasing), "Term"]
                }
                if(!is.null(pval)) enrichedRes[[i]] = enrichedRes[[i]][enrichedRes[[i]]$P.value < pval,]
                if(!is.null(padj)) enrichedRes[[i]] = enrichedRes[[i]][enrichedRes[[i]]$Adjusted.P.value < padj,]
                if(nrow(enrichedRes[[i]]) > 0 ) {
                        enrichedRes[[i]]$cluster = clusters[i]
                } else enrichedRes[[i]] =NULL
                Progress(i, length(clusters))
        }
        df_enriched <- data.table::rbindlist(enrichedRes) %>% as.data.frame()
        df_enriched = df_enriched[!is.na(df_enriched[, "Term"]),]
        df_enriched[," -log10(pval)"] = -log10(df_enriched$P.value)
        df_enriched[," -log10(padj)"] = -log10(df_enriched$Adjusted.P.value)
        if(verbose) print(round(nrow(df_enriched)/length(clusters)))
        
        if(isTRUE(Rowv) | isTRUE(Colv)) {
                mtx_enriched <- df_enriched[,c("Term","NES","cluster")]
                mtx_enriched %<>% tidyr::spread(cluster,NES)
                rownames(mtx_enriched) = mtx_enriched[,"Term"]
                mtx_enriched = mtx_enriched[,-grep("Term",colnames(mtx_enriched))]
                mtx_enriched %<>% as.matrix()
                mtx_enriched[is.na(mtx_enriched)] = 0
        }
        if(isTRUE(Rowv)) {
                hcr <- hclust(as.dist(1-cor(t(mtx_enriched), method="spearman")),
                              method="ward.D2")
                ddr <- as.dendrogram(hcr)
                rowInd <- order.dendrogram(ddr)
                order.row = rownames(mtx_enriched)[rowInd]
        } else {
                order.row = order.row[order.row %in% df_enriched[,"Term"]]
        }
        if(isTRUE(Colv)) {
                hcc <- hclust(as.dist(1-cor(mtx_enriched, method="spearman")),
                              method="ward.D2")
                ddc <- as.dendrogram(hcc)
                colInd <- order.dendrogram(ddc)
                order.col = colnames(mtx_enriched)[colInd]
        }
        df_enriched[,"Term"] %<>% as.factor
        df_enriched[,"Term"] %<>% factor(levels = order.row)
        # generate color pal_gsea scale based on NES range.
        rescale_colors <- function(colors = pal_gsea(), Range = range(df_enriched$NES)){
                if(Range[1]>0) return(pal_gsea()(12)[7:12])
                if(Range[2]<0) return(pal_gsea()(12)[1:6])
                if(Range[1]<0 & Range[2]>0) {
                        remove <- (Range[2] +Range[1]) / (Range[2] -Range[1])
                        if(remove>0) return(pal_gsea()(12)[max(1,12*remove):12])
                        if(remove<0) return(pal_gsea()(12)[1:(12+min(-1,12*remove)+1)])
                }
        }
        #font.ytickslab= min(font.ytickslab,round(height*300/dim(df_enriched)[1]))
        plot <- ggballoonplot(df_enriched, x = "cluster", y = "Term",
                              size = size, fill = fill,
                              size.range = c(1, 5),
                              font.ytickslab= font.ytickslab,
                              title = paste(pathway.name,"pathways enriched in",title),
                              legend.title = ifelse(fill =="NES",
                                                    "Normalized\nenrichment\nscore",
                                                    NULL),
                              xlab = "", ylab = "",...) +
                scale_fill_gradientn(colors = rescale_colors())+ #RPMG::SHOWPAL(ggsci::pal_gsea()(12))
                theme(plot.title = element_text(hjust = hjust,size = font.main))
        if(!is.null(order.col)) plot <- plot + scale_x_discrete(labels = order.col)
        if(size == "padj") plot = plot + scale_size(breaks=c(0,0.05,0.10,0.15,0.2,0.25),
                                                    labels=rev(c(0,0.05,0.10,0.15,0.2,0.25)))
        if(is.null(save_path)) {
                path <- paste0("output/",gsub("-","",Sys.Date()),"/")
                if(!dir.exists(path)) dir.create(path, recursive = T)
                save_path <- paste0(path,"Dotplot_",title,"_",pathway.name,
                                    "_",padj,"_",pval,".jpeg")
        }
        jpeg(save_path,units="in", width=width, height=height,res=600)
        print(plot)
        dev.off()
        if(do.return & return.raw) {
                return(enrichedRes)
        } else if(do.return) return(df_enriched)
}

df_enriched = EnrichDotPlot(stats=CSV_files, pathways=dbs[1:3], 
                        padj = 0.05,pval = 0.05,
                        order.by = c("3_BC","NES"),decreasing = F,
                        size = " -log10(padj)", fill = "Combined.Score",
                        pathway.name = "All curated gene sets",rotate.x.text = T,
                        title = "multiple lung cells",
                        font.xtickslab=12, font.main=17, font.ytickslab = 1,
                        font.legend = list(size = 15),font.label = list(size = 15),
                        do.return = T,
                        width = 8,height = 50)
write.csv(df_enriched, file = paste0("Yang/proximal_distal_terminal/Non-Integration/GSEA/cell types/Enrich_FDR0.25_pval0.05.csv"))
