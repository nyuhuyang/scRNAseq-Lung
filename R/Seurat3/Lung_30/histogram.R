# Libraries
invisible(lapply(c("Seurat","dplyr","magrittr","data.table","pbapply",
                   "tibble","tidyr","tidyverse","hrbrthemes","viridis","forcats"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
set.seed(101)

#======1.2 load  Seurat =========================
# load files
object = readRDS(file = "data/Lung_SCT_30_20200710.rds") 
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")

#X-axis: BC (all), IC (all), S, TASC, H, pre-C, C (all) - in this order
#Y-axis: log2(count+) gene expression in individual cells: KRT5, SERPINB3, MUC5AC, SCGB1A1, SCGB3A2, SFTPB, FOXA2, FOXJ1
#- establish max expression level for each of these genes across all samples and cells and consider it as 100%; then transform expression data as % from this 100% level.
#- then make the same type of plot you did, by using average transformed expression level (average % from max).
#- try to use smothering function for average lines (similar to how you did in pseudotime quantile graph)

#####################
#Colors: SERPINB3 - grey color or more distinct yellow 
#SFTPB the same red color as SCGB3A2

#Try KRT15 instead of KRT5 to get better shape
#Try CAPS instead of FOXJ1

gene_list <- list("Plot 1" = c("KRT15","SCGB1A1","SERPINB3","SFTPB"),
                  "Plot 2" = c("FOXA2","CAPS","MUC5AC","SCGB3A2"))

cbp3 <- list("Plot 1" = c("#0072B2", "#009E73", "#e69f00", "#ff0000"),
             "Plot 2" = c("#56B4E9","#000000","#329932","#ff0000"))
Idents(object) = "conditions"
conditions <- c("all_combined","distal","terminal","proximal","COPD")
for(i in seq_along(gene_list)){
        genes = gene_list[[i]]
        palette = cbp3[[names(gene_list)[i]]]
        
        for(con in conditions){
                if(con == "all_combined") {
                        sub_object <- object
                } else    sub_object <- subset(object, idents = con)
                Idents(sub_object) = "cell_types"
                exp = AverageExpression(sub_object,assays = "SCT",features = genes)
                exp = exp$SCT
                exp_max = apply(exp,1,max)
                
                exp = sweep(exp, 1, exp_max,"/")
                exp$gene = rownames(exp)
                data = gather(exp, "cell_types", "avg_logUMI", -gene)
                data$cell.type = plyr::mapvalues(data$cell_types,
                                                 from = c("BC","IC","S","TASC","H","p-C","C"),
                                                 to = 1:7)
                data$cell.type %<>% as.integer()
                
                g <- ggscatter(data = data,
                               x = "cell.type", y ="avg_logUMI",
                               color = "gene",
                               #shape = NA,                        # Extending the regression line
                               palette = palette,   
                               ylab = "max expression",
                               title = paste0(paste(gene_list[[i]], collapse = ", ")," in ",con),
                               add = "loess",add.params =  list(size=0.5),
                               conf.int.level = 0.15,
                               conf.int = T)+TitleCenter()+
                        scale_x_continuous(name ="cell_types",
                                           breaks = 1:7,
                                           labels = c("BC","IC","S","TASC","H","p-C","C"))
                jpeg(paste0(path,names(gene_list)[i],"_",con,".jpeg"), units="in", width=10, height=7,res=600)
                print(g)
                dev.off()
                Progress(which(conditions %in% con), length(conditions))
        }
}



cell_types <- c("BC1","BC2","IC1","IC2","IC3","S","TASC","H","p-C","C1","C2","C3")
(genes <- sort(c("KRT5", "SERPINB3", "MUC5AC", "SCGB1A1", "SCGB3A2", "SFTPB", "FOXA2", "FOXJ1")))

Idents(object) = "cell_types"
object <- subset(object, idents = cell_types)
object$cell_types %<>% gsub("1|2|3","",.)
object$cell_types %<>% factor(levels = c("BC","IC","S","TASC","H","p-C","C"))

Idents(object) = "cell_types"
exp = AverageExpression(object,assays = "SCT",features = genes)
exp = exp$SCT
exp_max = apply(exp,1,max)

exp = sweep(exp, 1, exp_max,"/")
exp$gene = rownames(exp)
data = gather(exp, "cell_types", "avg_logUMI", -gene)
data$cell.type = plyr::mapvalues(data$cell_types,
                                  from = c("BC","IC","S","TASC","H","p-C","C"),
                                  to = 1:7)
data$cell.type %<>% as.integer()
#=========
#Exp = FetchData(object,vars = c(genes,"cell_types","conditions"))

#data1 <- Exp %>% gather("gene", "avg_logUMI", -c("cell_types","conditions")) %>%
#        group_by(gene) %>% 
#        mutate(max_exp =max(avg_logUMI)) %>%
#        mutate(scaled_exp = avg_logUMI/max_exp) %>%
#        ungroup()
#data1$cell.type = plyr::mapvalues(data1$cell_types,
#                                 from = c("BC","IC","S","TASC","H","p-C","C"),
#                                  to = 1:7)
#data1$cell.type %<>% as.integer()
#===
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
levels = c(0.15,0.5,0.95)
for(lvl in levels){
        g <- ggscatter(data = data,
                  x = "cell.type", y ="avg_logUMI",
                  color = "gene",
                  #shape = NA,                        # Extending the regression line
                  palette = cbp1,  
                  ylab = "max expression",
                  title = paste0("all_combined_conf.int.level=",lvl),
                  add = "loess",add.params =  list(size=0.5),
                  conf.int.level = lvl,
                  conf.int = T)+TitleCenter()+
                scale_x_continuous(name ="cell_types",
                                   breaks = 1:7,
                                   labels = c("BC","IC","S","TASC","H","p-C","C"))
        jpeg(paste0(path,"gene_curves_all_combined_conf.int.level=",lvl,".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
}




#=========== for each Region =========
Idents(object) = "conditions"
conditions <- c("distal","terminal","proximal","COPD")
for(con in conditions){
        sub_object <- subset(object, idents = con)
        Idents(sub_object) = "cell_types"
        exp = AverageExpression(sub_object,assays = "SCT",features = genes)
        exp = exp$SCT
        exp_max = apply(exp,1,max)
        
        exp = sweep(exp, 1, exp_max,"/")
        exp$gene = rownames(exp)
        data = gather(exp, "cell_types", "avg_logUMI", -gene)
        data$cell.type = plyr::mapvalues(data$cell_types,
                                         from = c("BC","IC","S","TASC","H","p-C","C"),
                                         to = 1:7)
        data$cell.type %<>% as.integer()

        g <- ggscatter(data = data,
                  x = "cell.type", y ="avg_logUMI",
                  color = "gene",
                  #shape = NA,                        # Extending the regression line
                  palette = cbp1,   
                  ylab = "max expression",
                  title = paste0(con,"_conf.int.level=0.15"),
                  add = "loess",add.params =  list(size=0.5),
                  conf.int.level = 0.15,
                  conf.int = T)+TitleCenter()+
                scale_x_continuous(name ="cell_types",
                                   breaks = 1:7,
                                   labels = c("BC","IC","S","TASC","H","p-C","C"))
        jpeg(paste0(path,"gene_curves_",con,"_conf.int.level=0.15.jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        Progress(which(conditions %in% con), length(conditions))
}

#####===================================================

object@meta.data$cell_types %<>% factor(levels = cell.type)
hg_19_38 = readRDS("data/RNA-seq/hg_19_38.rds")
table(rownames(object) %in% hg_19_38$hg38)
table(rownames(object) %in% c(hg_19_38$hg19,hg_19_38$hg38))


cell.type_list <- list("Epithelial" = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","TASC",
                                        "H","p-C","C1","C2","C3","Ion","NE","ME","g-Muc",
                                        "g-Ser","AT1","AT2","AT2-1","AT2-p"),
                       "Stromal"=c("F1","F2","F3","F4","Cr","Gli","Nr","SM1",
                                   "SM2","SM3","Pr","En-a","En-c","En-c1","En-v","En-l","En-sm","En-p"),
                       "Immune" = c("MC","Neu","Mon","M0","M1","M1-2","M2","M-p","DC","p-DC",
                                    "B","PC","T-cn","T-reg","T-int","T-rm","T-NK","T-ifn","T-p"))
cell.type <- unlist(cell.type_list,use.names = FALSE)
superfamily <- c("Epithelial","Stromal","Immune")

feature_list <- list(
        "Panel_1" = c("IFNG", "IFNGR1", "IFNGR2", "TNF", "TNFRSF1A"),
        "Panel_2" = c("IL33", "IL1RL1", "TSLP", "CRLF2", "IL3", "IL5", "IL13", "IL25"),
        "Panel_3" = c("AREG", "KIT", "EREG", "EGFR", "EGFR", "ERBB2", "ERBB4", "NMU", "NMUR1"),
        "Panel_4" = c("GATA3", "RORC", "CRLF2", "PTGDR2", "KLRG1", "BCL11B", "IL7R", "IL7")
)
features <- unlist(feature_list,use.names = FALSE) %>% unique
table(features %in% rownames(object))
features[!(features %in% rownames(object))]
# create expression data frame
exp = FetchData(object,vars = c(unique(features),"cell_types"))
exp %<>% gather(key="gene", value="expression", -cell_types)

exp$cell_types %<>% factor(levels = cell.type)
sub_exp <- exp %>% filter(gene %in% feature_list[["Panel_1"]]) %>%
        filter(cell_types %in% cell.type_list[["Epithelial"]])

g <- sub_exp %>% filter(expression > 0) %>%
        ggplot( aes(x=expression, color=gene, fill=gene)) +
        geom_density(alpha=0.6, binwidth = 0.1) +
        scale_fill_viridis(discrete=TRUE) +
        scale_color_viridis(discrete=TRUE) +
        theme_ipsum() +
        theme(
                legend.position="none",
                panel.spacing = unit(0.1, "lines"),
                axis.text.x = element_text(angle = 0,size = 10,
                                           vjust = 0.5, hjust=1)
                #panel.grid.major = element_blank(),
                #panel.grid.minor = element_blank()
        ) +
        xlab("") +
        ylab("gene expression counts") +
        facet_grid(gene~cell_types)
jpeg(paste0(path, "Panel_1","_","Epithelial_density>0",".jpeg"), units="in", width=15, height=7,res=600)
print(g)
dev.off()

object <- subset(object, idents = cell.type)
TASC <- subset(object, idents = "TASC")

object@meta.data$cell_types %<>% factor(levels = cell.type)

jpeg(paste0(path, "Allcells_SCGB1A1.jpeg"), units="in", width=7, height=10,res=600)
RidgePlot(object, features = "SCGB1A1",group.by = "cell_types", sort = TRUE,ncol = 1) + NoLegend()
dev.off()

TASC@meta.data$conditions %<>% factor(levels = rev(c("proximal","distal","terminal","COPD")))
jpeg(paste0(path, "TASC_SCGB1A1_SFTPA1_SFTPA2.jpeg"), units="in", width=10, height=7,res=600)
RidgePlot(TASC, features = c("SCGB1A1","SFTPA1","SFTPA2"),ncol = 2,
          cols = rev(c("#1F78B4","#4ca64c","#E6AB02",'#FF4136')),
          group.by = "conditions") + NoLegend()
dev.off()

TASC@meta.data$conditions %<>% factor(levels = c("proximal","distal","terminal","COPD"))
Cols =c("#1F78B4","#4ca64c","#E6AB02",'#FF4136')

plot1 <- FeatureScatter(TASC, feature1 = "SCGB1A1", feature2 = "SFTPA1",group.by = "conditions",cols=Cols)
plot2 <- FeatureScatter(TASC, feature1 = "SCGB1A1", feature2 = "SFTPA2",group.by = "conditions",cols=Cols)
plot3 <- FeatureScatter(TASC, feature1 = "SCGB1A1", feature2 = "SCGB3A2",group.by = "conditions",cols=Cols)
plot4 <- FeatureScatter(TASC, feature1 = "SCGB1A1", feature2 = "SFTPB",group.by = "conditions",cols=Cols)

jpeg(paste0(path, "TASC_FeatureScatter.jpeg"), units="in", width=10, height=7,res=600)
((plot1 & NoLegend())+ plot2 + (plot1 & NoLegend())+ plot4) 
dev.off()

Features = c("SFTPA1","SFTPA2","SCGB3A2","SFTPB")
exp = FetchData(TASC,vars = c("SCGB1A1", Features))
UMI = colMeans(exp)
cor_rest = sapply(Features,function(x){
       cor.test(exp$SCGB1A1, exp[,x],method = c("pearson"))$estimate[["cor"]]
})

p_rest = sapply(Features,function(x){
        cor.test(exp$SCGB1A1, exp[,x], method = c("pearson"))$p.value
})

df = bind_rows(list(UMI, cor_rest, p_rest)) %>% as.data.frame()

rownames(df) = c("avg_UMI","corrlation","p value")
write.csv(df,file = paste0(path,"SCGB1A1_correlation.csv"))

# histogram
#X-axis: BC (all), IC (all), S, TASC, H, pre-C, C (all) - in this order
#Y-axis: log2(count+) gene expression in individual cells: KRT5, SERPINB3, MUC5AC, SCGB1A1, SCGB3A2, SFTPB, FOXA2, FOXJ1
cell_types <- c("BC1","BC2","IC1","IC2","IC3","S","TASC","H","p-C","C1","C2","C3")
genes <- c("KRT5", "SERPINB3", "MUC5AC", "SCGB1A1", "SCGB3A2", "SFTPB", "FOXA2", "FOXJ1")

Idents(object) = "cell_types"
object <- subset(object, idents = cell_types)
object$cell_types %<>% gsub("1|2|3","",.)
object$cell_types %<>% factor(levels = c("BC","IC","S","TASC","H","p-C","C"))
Idents(object) = "cell_types"

RidgePlot(object, features = genes[1],ncol = 1,group.by = "cell_types"
          #cols = #rev(c("#1F78B4","#4ca64c","#E6AB02",'#FF4136'))
          ) + NoLegend()

jpeg(paste0(path, "Dot_plot_all_combined.jpeg"), units="in", width=4.5, height=2.9,res=600)
DotPlot.1(objectv,features = rev(genes),exp.max = 7,dot.scale = 6, 
          scale = FALSE,cluster.idents = FALSE,cluster.features = FALSE,
          cols = c("blue","green","yellow","orange","chocolate1","red"))+
        ggtitle("all groups combined")+
        theme(axis.line=element_blank(),
              text = element_text(size=14),#16
              panel.grid.major = element_blank(),
              axis.text.x = element_text(size=14,#16
                                         angle = 90, 
                                         hjust = 0,
                                         vjust= 0.5),
              legend.title = element_text(size=12), #change legend title font size
              legend.text = element_text(size=12),
              legend.key.size = unit(0.25, 'cm'), #change legend key size
              legend.key.width = unit(0.25, 'cm'), #change legend key width
              axis.title.x = element_blank(),
              axis.text.y = element_text(size=14),#14
              axis.title.y = element_blank(),
              plot.title = element_text(face = "plain",hjust = 0.5))+
        scale_y_discrete(position = "right")+
        coord_flip()
dev.off()

Idents(object) = "conditions"
conditions <- c("distal","terminal","proximal","COPD")
for(con in conditions){
        sub_object <- subset(object, idents = con)
        Idents(sub_object) = "cell_types"
        jpeg(paste0(path, "Dot_plot_",con,".jpeg"), units="in", width=4.5, height=2.9,res=600)
        g <- DotPlot.1(sub_object,features = rev(genes),exp.max = 7,dot.scale = 6, 
                  scale = FALSE,cluster.idents = FALSE,cluster.features = FALSE,
                  cols = c("blue","green","yellow","orange","chocolate1","red"))+
                ggtitle(con)+
                theme(axis.line=element_blank(),
                      text = element_text(size=14),#16
                      panel.grid.major = element_blank(),
                      axis.text.x = element_text(size=14,#16
                                                 angle = 90, 
                                                 hjust = 0,
                                                 vjust= 0.5),
                      legend.title = element_text(size=12), #change legend title font size
                      legend.text = element_text(size=12),
                      legend.key.size = unit(0.25, 'cm'), #change legend key size
                      legend.key.width = unit(0.25, 'cm'), #change legend key width
                      axis.title.x = element_blank(),
                      axis.text.y = element_text(size=14),#14
                      axis.title.y = element_blank(),
                      plot.title = element_text(face = "plain",hjust = 0.5))+
                scale_y_discrete(position = "right")+
                coord_flip()
        print(g)
        dev.off()
        Progress(which(conditions %in% con), length(conditions))
}


for(gene in genes){
        jpeg(paste0(path, "VlnPlot/","all_combined_",gene,".jpeg"), units="in", width=10, height=7,res=600)
        print(VlnPlot(object = object, features = gene))
        dev.off()
}

Idents(object) = "conditions"
conditions <- c("distal","terminal","proximal","COPD")
for(con in conditions){
        sub_object <- subset(object, idents = con)
        Idents(sub_object) = "cell_types"
        for(gene in genes){
                jpeg(paste0(path, "VlnPlot/",con,"_",gene,".jpeg"), units="in", width=10, height=7,res=600)
                print(VlnPlot(object = sub_object, features = gene))
                dev.off()
        }
        Progress(which(conditions %in% con), length(conditions))
}

