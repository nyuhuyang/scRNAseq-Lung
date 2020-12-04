########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","kableExtra",
                   "magrittr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

set.seed(101)
#======1.2 load  Seurat =========================
# load files
object = readRDS(file = "data/Lung_30_20200710.rds") 
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
Idents(object) = "conditions"

df_abbr <- readxl::read_excel("doc/Annotations/Cell type abbreviation.xlsx",
                              sheet = "abbreviation_20201116")

table(df_abbr$Abbreviation %in% object$annotations3)
object$long_annotations3 <- plyr::mapvalues(object$annotations3,
                                            from = df_abbr$Abbreviation,
                                            to = df_abbr$`Cell types`)
# distal,terminal,proximal
sub_object <- subset(object, idents = c("distal","terminal","proximal"))
Idents(sub_object) = "conditions"
Idents(sub_object) %<>% factor(levels = c("proximal", "distal", "terminal"))
Idents(sub_object) %<>% factor(levels = c("terminal","distal", "proximal"))

UMAPPlot.1(sub_object,group.by = "conditions", cols = c("#4ca64c","#1F78B4","#E6AB02"), 
           label = F, label.repel = F, 
           pt.size = 1,alpha.by = "conditions",
           width=10, height=10, 
           no.legend = F,do.return = F, do.print = T)

condition_annotaion <- table(sub_object$long_annotations3, sub_object$conditions) %>% 
    as.data.frame.matrix()
condition_annotaion %>% kable() %>% kable_styling()
condition_annotaion %<>% sweep(2, colSums(.),"/") %>%
    sweep(1, rowSums(.),"/")

Rowv = F
if(isTRUE(Rowv)){
    hcr <- hclust(as.dist(1-cor(t(condition_annotaion), method="spearman")),
                  method="ward.D2")
    ddr <- as.dendrogram(hcr)
    rowInd <- order.dendrogram(ddr)
    condition_annotaion = condition_annotaion[rowInd,]
} else condition_annotaion %<>% .[order(.$proximal,decreasing = T),]


jpeg(paste0(path,"region_specific_cell_types.jpeg"), units="in", width=5, height=10,res=600)
ggballoonplot(condition_annotaion, 
              fill = "value",
              size.range = c(1, 5),
              rotate.x.text = FALSE,
              font.xtickslab = 12)+
    scale_fill_viridis_c(option = "C")
dev.off()

# distal,terminal,proximal, COPD
UMAPPlot.1(object,group.by = "conditions", cols = c("#F0027F","#4ca64c","#1F78B4","#E6AB02"), 
           label = F, label.repel = F, 
           pt.size = 1,alpha.by = "conditions",
           width=10, height=10, 
           no.legend = T,do.return = F, do.print = T)

condition_annotaion <- table(object$long_annotations3, object$conditions) %>% 
    as.data.frame.matrix()
condition_annotaion %>% kable() %>% kable_styling()
condition_annotaion %<>% sweep(2, colSums(.),"/") %>%
    sweep(1, rowSums(.),"/")

condition_annotaion %<>% .[order(.$proximal,decreasing = T),]
condition_annotaion %<>% .[,c("proximal", "distal", "COPD", "terminal")]
jpeg(paste0(path,"region_specific_cell_types~.jpeg"), units="in", width=6, height=10,res=600)
ggballoonplot(condition_annotaion, 
              fill = "value",
              size.range = c(1, 5),
              rotate.x.text = FALSE,
              font.xtickslab = 12)+
    scale_fill_viridis_c(option = "C")
dev.off()


Idents(object) = "long_annotations3"
s_d <- subset(object, idents = "Secretory cells-Distal")
s_d@meta.data$conditions %<>% as.factor()
s_d@meta.data$conditions %<>% factor(levels = c("proximal", "distal", "COPD", "terminal"))
DotPlot(s_d, features = c("SCGB1A1","SCGB3A1","SCGB3A2","SFTPB","MUC5B"),
        split.by = "conditions",cols = c("#4ca64c","#E6AB02","#1F78B4","#F0027F"))+ 
    coord_flip() + RotatedAxis()

# Color of dots (blue, green, yellow, red; see example of color range in Figure 5a in the attached paper) - expression level.
path = "Yang/Lung_30/Dotplots/No Clustering/"
Idents(object) = "annotations3"

groups = c("Epithelial","Structural","Immune")
for(i in 1:length(groups)){
    df <- readxl::read_excel("doc/30 DEGs for dot plots - revised 12-1-20.xlsx",
                             sheet = groups[i])
    cell.types = colnames(df)[2:ncol(df)]
    genes = df$...1
    sub_object <- subset(object, idents = cell.types)
    sub_object@meta.data$annotations3 %<>% as.factor()
    sub_object@meta.data$annotations3 %<>% factor(levels = cell.types)
    Idents(sub_object) = "annotations3"
    
    options <- c("Max_6","Max_7","Max_auto")[2]
    for(k in 1:length(options)){
        save.path <- paste0(path,options[k],"/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        exp.max = switch(options[k],
                         "Max_6" = 6,
                         "Max_7" = 7,
                         "Max_auto" = NA)
        g <- DotPlot.1(sub_object, features = rev(genes),
                       log.data = log2,exp.max = exp.max,
                       scale = TRUE,
                       cluster.idents = FALSE,
                       cluster.features = FALSE,
                       cols = c("blue","green","yellow","orange","chocolate1","red"))+
            ggtitle(groups[i])+ TitleCenter()+
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                  text = element_text(size=16),
                  axis.text.x = element_text(size=16,
                                             angle = 90, 
                                             hjust = 1),
                  axis.text.y = element_text(size=16),
                  
                  plot.title = element_text(face = "plain"))+
            coord_flip()
        width = switch(groups[i],
                       "Epithelial" = 9,
                       "Structural" = 7.5,
                       "Immune" = 8.5)
        jpeg(paste0(save.path,"Dotplot_",groups[i],".jpeg"), family = "Arial",
             units="in", width=width, height=10,res=900)
        print(g)
        dev.off()
    }
}

