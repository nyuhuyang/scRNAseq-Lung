########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","tidyr","magrittr","ggpubr","MAST",
                   "gplots","fgsea"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#######################
# box plot
#######################

# for major cell ===========================
# load data
cell_dist <- readxl::read_excel("doc/Cell distribution - for plots.xlsx")
cell_dist = cell_dist[!is.na(cell_dist$Proximal),]
colnames(cell_dist) = gsub(" ",".",cell_dist[1,])
cell_dist = cell_dist[-1,]
cell_dist = cell_dist[,-grep("Families|Code",colnames(cell_dist))]

regions = c("Proximal","Distal","Terminal")
major_cells <- c("Epithelial","Stromal","Endothelial","Immune","Neuronal")
major_cell_dist = as.data.frame(cell_dist[1:5,])
df_major_cell_dist <- gather(major_cell_dist, key = "samples", value = "percentage",-Cell.categories)
head(df_major_cell_dist)
df_major_cell_dist$percentage %<>% as.numeric() 
df_major_cell_dist$regions = gsub("-R","",df_major_cell_dist$samples)
df_major_cell_dist$regions %<>% gsub(".*-","",.)
df_major_cell_dist$regions %<>% plyr::mapvalues(from = c("P","D","T"),
                                       to = regions)
df_major_cell_dist$regions %<>% as.factor()
df_major_cell_dist$regions %<>% factor(levels = rev(regions))
table(df_major_cell_dist$regions)
#df_major_cell_dist = df_major_cell_dist[,-grep("samples",colnames(df_major_cell_dist))]
#rownames(df_major_cell_dist) = df_major_cell_dist$Cell.categories
df_major_cell_dist1 <- spread(df_major_cell_dist, Cell.categories,percentage, fill = 0)
head(df_major_cell_dist1)

jpeg(paste0(path,"Cell distribution.jpeg"), units="in", width=10, height=7,res=600)
ggboxplot(df_major_cell_dist1, x = "regions",
          y = rev(major_cells),
          merge = "flip",
          ylab = "Cell percenrage", 
          xlab = "",
          palette = c('#3D9970','#FF4136','#FF851B'),
          color = "regions",
          add = "jitter",
          #add.params = list(size = 0.1, jitter = 0.2),
          rotate = TRUE)
dev.off()

my_comparisons <- list(c("Proximal", "Distal"), c("Distal", "Terminal"),c("Proximal", "Terminal"))

for(c in major_cells){
        g <- ggboxplot(df_major_cell_dist1, x = "regions",
                       y = c,
                       merge = "flip",
                       ylab = "Cell percenrage", 
                       xlab = "",
                       palette = c('#3D9970','#FF4136','#FF851B'),
                       color = "regions",
                       add = "jitter",
                       #add.params = list(size = 0.1, jitter = 0.2),
                       rotate = TRUE)+
                stat_compare_means(comparisons = my_comparisons)
        jpeg(paste0(path,"Cell distribution-",c,".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
}

# for minor cell =======================
# load data
cell_dist <- readxl::read_excel("doc/Cell distribution - for plots.xlsx")
cell_dist = cell_dist[!is.na(cell_dist$Proximal),]
colnames(cell_dist) = gsub(" ",".",cell_dist[1,])
cell_dist = cell_dist[-1,]
cell_dist = cell_dist[,-grep("Cell.categories",colnames(cell_dist))]
#https://stackoverflow.com/questions/10554741/fill-in-data-frame-with-values-from-rows-above/32536507
cell_dist$Families = zoo::na.locf(cell_dist$Families)
cell_dist = as.data.frame(cell_dist[-c(1:5),])

regions = c("Proximal","Distal","Terminal")
(minor_cells <- unique(cell_dist$Families))
for(cell in minor_cells) {
        print(cell)
        cell_path <- paste0(path,cell,"/")
        if(!dir.exists(cell_path))dir.create(cell_path, recursive = T)
        
        minor_cell_dist = cell_dist[cell_dist$Families %in% cell,-1]
        df_minor_cell_dist <- gather(minor_cell_dist, key = "samples", value = "percentage",-Code)
        head(df_minor_cell_dist,2)
        df_minor_cell_dist$percentage %<>% as.numeric() 
        df_minor_cell_dist$regions = gsub("-R","",df_minor_cell_dist$samples)
        df_minor_cell_dist$regions %<>% gsub(".*-","",.)
        df_minor_cell_dist$regions %<>% plyr::mapvalues(from = c("P","D","T"),
                                                        to = regions)
        df_minor_cell_dist$regions %<>% as.factor()
        df_minor_cell_dist$regions %<>% factor(levels = rev(regions))
        table(df_minor_cell_dist$regions)
        df_minor_cell_dist1 <- spread(df_minor_cell_dist, Code,percentage, fill = 0)
        ggboxplot(df_minor_cell_dist1, x = "regions",
                  y = rev(minor_cell_dist$Code),
                  bxp.errorbar = TRUE,
                  merge = "flip",
                  width = 0.7,
                  size = NULL,
                  ylab = "Cell percenrage", 
                  xlab = "",
                  palette = c('#3D9970','#FF4136','#FF851B'),
                  color = "regions",
                  add = "jitter",
                  #add.params = list(size = 0, jitter = 0),
                  rotate = TRUE)
        jpeg(paste0(cell_path,"Cell distribution-",cell,".jpeg"), units="in", width=10, height=20,res=600)
        
        print(g)
        dev.off()
        
        my_comparisons <- list(c("Proximal", "Distal"), c("Distal", "Terminal"),c("Proximal", "Terminal"))
        colnames(df_minor_cell_dist1) %<>% gsub("-|/","_",.)
        (minor_cell_types <- colnames(df_minor_cell_dist1)[-c(1:2)])

        for(c in minor_cell_types){
                df_minor_cell_dist2 <- df_minor_cell_dist1[,c("samples", "regions",c)]
                g1 <- ggboxplot(df_minor_cell_dist1, x = "regions",
                               y = c,
                               merge = "flip",
                               ylab = "Cell percenrage", 
                               xlab = "",
                               palette = c('#3D9970','#FF4136','#FF851B'),
                               color = "regions",
                               add = "jitter",
                               #add.params = list(size = 0.1, jitter = 0.2),
                               rotate = TRUE)+
                        stat_compare_means(comparisons = my_comparisons)
                jpeg(paste0(cell_path,"Cell distribution-",c,".jpeg"), units="in", width=10, height=7,res=600)
                print(g1)
                dev.off()
        }
        Progress(which(minor_cells %in% cell),length(minor_cells))
}

#######################
# GSEA
#######################
Yang_path <- "Yang/proximal_distal_terminal/Non-Integration/DEGs/cell_types/"
csv_files <- list.files(Yang_path,pattern = "Lung_24-FC0_markers_") 
csv_num <- gsub("Lung_24-FC0_markers_","",csv_files) %>% 
        gsub("_.*","",.) %>%  as.numeric()
num <- 1:67
num[!(num %in% csv_num)]
genes.de <- list()
for (i in 1:length(x = csv_files)) {
        genes.de[[i]] <- read.csv(paste0(Yang_path,csv_files[i]), )
}
gde.all <- data.frame()
idents.all <- gsub(".*_","",csv_files) %>% gsub(".csv","",.)
test.use = "MAST";node = NULL; return.thresh = 1
for (i in 1:length(x = idents.all)) {
        if (is.null(x = unlist(x = genes.de[i]))) {
                next
        }
        gde <- genes.de[[i]]
        if (nrow(x = gde) > 0) {
                if (is.null(x = node) || test.use %in% c("bimod", "t")) {
                        gde <- gde[order(gde$p_val, -gde[, 2]), ]
                        gde <- subset(x = gde, subset = p_val < return.thresh)
                }
                if (nrow(x = gde) > 0) {
                        gde$cluster <- idents.all[i]
                        gde$gene <- gde$X
                }
                if (nrow(x = gde) > 0) {
                        gde.all <- rbind(gde.all, gde)
                }
        }
        Progress(i,length(x = idents.all))
}
rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
write.csv(gde.all,file = paste0(Yang_path,"Lung_24-FC0_cell_types.csv"))
# get Pathways
hallmark <- gmtPathways("../seurat_resources/msigdb/h.all.v7.0.symbols.gmt")
kegg <- gmtPathways("../seurat_resources/msigdb/c2.cp.kegg.v7.0.symbols.gmt")
bp <- gmtPathways("../seurat_resources/msigdb/c5.bp.v7.0.symbols.gmt")
mf <- gmtPathways("../seurat_resources/msigdb/c5.mf.v7.0.symbols.gmt")
go <- c(bp, mf)
# Now, run the fgsea algorithm with 1000 permutations:
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(kegg) = gsub("KEGG_","",names(kegg))
names(go) = gsub("GO_","",names(go))
pathways_list <- list("HALLMARK"=hallmark,"KEGG"= kegg, "GO"= go)

fgseaRes = FgseaDotPlot(stats=gde.all, pathways=pathways_list[[1]], nperm=1000,
                        padj = 0.25,pval = 0.05,
                        Rowv = T,Colv = T,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = names(pathways_list[1]),rotate.x.text = T,
                        title = "multiple lung cells",
                        font.xtickslab=10, font.main=17, font.ytickslab = 11,
                        font.legend = list(size = 16),font.label = list(size = 16),
                        do.return = T,
                        width = 15,height = 10)
write.csv(fgseaRes, file = paste0(path,"Hallmark_FDR0.25_pval0.05.csv"))

fgseaRes = FgseaDotPlot(stats=gde.all, pathways=pathways_list[[2]], nperm=1000,
                        padj = 0.25,pval = 0.05,
                        Rowv = T,Colv = T,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = names(pathways_list[2]),rotate.x.text = T,
                        title = "multiple lung cells",
                        font.xtickslab=10, font.main=17, font.ytickslab = 7,
                        font.legend = list(size = 14),font.label = list(size = 14),
                        do.return = T,
                        width = 15,height = 15)
write.csv(fgseaRes, file = paste0(path,"KEGG_FDR0.25_pval0.05.csv"))

fgseaRes = FgseaDotPlot(stats=gde.all, pathways=pathways_list[[3]], nperm=1000,
                        padj = 0.25,pval = 0.05,
                        Rowv = T,Colv = T,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = names(pathways_list[3]),rotate.x.text = T,
                        title = "multiple lung cells",
                        font.xtickslab=10, font.main=17, font.ytickslab = 1,
                        font.legend = list(size = 14),font.label = list(size = 14),
                        do.return = T,
                        width = 15,height = 30)
write.csv(fgseaRes, file = paste0(path,"GO_FDR0.25_pval0.05.csv"))

#######################
# Dendrogram
#######################
(load(file = paste0("data/Lung_24_20191206.Rda")))
DefaultAssay(object)  = "SCT"
Idents(object) = "cell.types"
object %<>% sortIdent
exp <- AverageExpression(object, assays = "SCT")
top_genes <- head(VariableFeatures(object), 2000)
exp =  exp$SCT[top_genes,]
exp = exp %>% t %>% scale %>% t
save(exp, file = paste0(path,"exp.Rda"))
load(file = paste0(path,"exp.Rda"))
write.csv(exp,paste0(path,"AverageExpression.csv"))

write.csv(exp,paste0(path,"AverageExpression_zscore.csv"))
methods_list <- list(c("spearman","complete",T),
                     c("spearman","complete",F),
                     c("spearman","ward.D2",T),
                     c("spearman","ward.D2",F))
#                     c("spearman","average"),
#                     c("pearson","complete"),
#                     c("pearson","ward.D2"),
#                     c("pearson","average"))
for(i in seq_along(methods_list)) {
        heatmap_path <- paste0(path,"Heatmap2_",methods_list[[i]][1],
                               "_",methods_list[[i]][2],
                               "_RowMean=",methods_list[[i]][3],"/")
        if(!dir.exists(heatmap_path)) dir.create(heatmap_path, recursive = T)
        # column hcluster
        hc <- hclust(as.dist(1-cor(exp, method=methods_list[[i]][1])),
                     method=methods_list[[i]][2])
        cc = as.character(as.numeric(as.factor(hc$labels)))
        # row hcluster
        hc_genes <- hclust(as.dist(1-cor(t(exp), method=methods_list[[i]][1])),
                     method=methods_list[[i]][2])
        jpeg(paste0(heatmap_path,"Heatmap2_",methods_list[[i]][1],
                    "_",methods_list[[i]][2],".jpeg"),
             units="in", width=10, height=7,res=600)
        g <- heatmap.2(exp,
                       Colv = as.dendrogram(hc),
                       Rowv = if(methods_list[[i]][3]) {T} else {as.dendrogram(hc_genes)},
                       ColSideColors = cc, 
                       trace ="none",
                       dendrogram = "both",
                       key.xlab = "scale log nUMI",
                       cexRow = 0.5,
                       margins = c(13,5),
                       breaks = seq(-3,3,length.out = 101),
                       col = bluered,
                       main = paste("Clustering dendrogram for all cell types\n based on top 2000", 
                                    "variable genes\n by",methods_list[[i]][1],
                                    "correlation",methods_list[[i]][2], "linkage"))
        dev.off()
        write.csv(rownames(exp)[rev(g$rowInd)],
                  file = paste0(heatmap_path,"Heatmap2_celltype_",methods_list[[i]][1],
                                "_",methods_list[[i]][2],"_rowInd.csv"))
        write.csv(colnames(exp)[g$colInd],
                  file = paste0(heatmap_path,"Heatmap2_celltype_",methods_list[[i]][1],
                                "_",methods_list[[i]][2],"_colInd.csv"))
        write.csv(capture.output(str(as.dendrogram(hc_genes))),
                  file = paste0(heatmap_path,"Tree_genes_",methods_list[[i]][1],
                                "_",methods_list[[i]][2],".csv"))
        write.csv(capture.output(str(as.dendrogram(hc))),
                  file = paste0(heatmap_path,"Tree_celltypes_",methods_list[[i]][1],
                                "_",methods_list[[i]][2],".csv"))
        group.use = gsub(":.*","",colnames(exp)[g$colInd])
        group.use = factor(group.use,levels = unique(group.use))
        
        DoHeatmap.matrix(object = exp, features = rownames(exp)[rev(g$rowInd)],
                         cells = colnames(exp)[g$colInd],
                         group.by = group.use,
                         group.colors = Singler.colors,
                         size = 4, cex.row=1, draw.lines = F,raster = F,
                         title = paste(methods_list[[i]][1],methods_list[[i]][2]),
                         save.path = heatmap_path,
                         pal_gsea = F, do.print = T, no.legend = T, 
                         angle = 90,width = 10,height=30)
        Progress(i, length(methods_list))
}

