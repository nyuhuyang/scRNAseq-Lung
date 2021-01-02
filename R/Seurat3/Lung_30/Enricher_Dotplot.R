########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("dplyr","magrittr","data.table","pbapply","tibble","tidyr","patchwork"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

set.seed(101)
# ========== define dotplot rownames and colnames =================
cell.type_list <- list("Epithelial" = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","d-S",
                                        "H","p-C","C1","C2","C3","Ion","NE","ME","g-Muc",
                                        "g-Ser","AT1","AT2","AT2-1","AT2-p"),
                       "Structural"=c("F1","F2","F3","F4","Cr","Gli","Nr","SM1",
                                      "SM2","SM3","Pr","En-a","En-c","En-c1","En-v","En-l","En-sm","En-p"),
                       "Immune" = c("MC","Neu","Mon","M0","M1","M1-2","M2","M-p","DC","p-DC",
                                    "B","PC","T-cn","T-reg","T-int","T-rm","T-NK","T7","T-p"))
selected_tissues <- list("GTEx tissue" = c("lung",
                                           "esophagus",
                                           "skin",
                                           "stomach",
                                           "salivary gland",
                                           "testis",
                                           "brain",
                                           "adipose tissue",
                                           "blood vessel",
                                           "nerve",
                                           "blood",
                                           "spleen"),
                         "Human Gene Atlas" = c("Trachea",
                                                "BronchialEpithelialCells",
                                                "Tongue",
                                                "OlfactoryBulb",
                                                "CD33+ Myeloid",
                                                "CD14+ Monocytes",
                                                "BDCA4+ DendriticCells",
                                                "CD19+ BCells(neg. sel.)",
                                                "CD4+ Tcells",
                                                "CD8+ Tcells",
                                                "CD56+ NKCells"))
csv_list <- c("enrichR_EVG_GTEx_Tissue_Sample_Gene_Expression_Profiles_up.csv",
              "enrichR_EVG_Human_Gene_Atlas.csv")
adj = c(10^-4,10^-5)[1]
save.path = paste0(path,"adj < ",adj[1],"/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

superfamily <- c("Epithelial","Structural","Immune")
#======1.2 load  Seurat =========================
# load files
object = readRDS(file = "data/Lung_30_20200710.rds") 
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")

anno <- readxl::read_excel("doc/Annotations/Cell type abbreviation.xlsx")
table(anno$Abbreviation %in% object$annotations3)
object$cell_types <- plyr::mapvalues(object$annotations3,
                                     from = anno$Abbreviation,
                                     to = anno$`Revised abbreviations`)
Idents(object) = "cell_types"
dotplot_df <- readxl::read_excel("doc/Instructions - expanded dotplot - cell types_Yang modified.xlsx")
features_list <- dotplot_df[1:40,superfamily]

# for top 3 figures ===========
g <- list()
for(i in seq_along(superfamily)){
    group = superfamily[i]
    sub_object <- subset(object, idents =  cell.type_list[[group]])
    Idents(sub_object) %<>% factor(levels = cell.type_list[[group]])
    features = features_list[[group]]
    g[[i]] <- DotPlot.1(sub_object, features = rev(features),
                   log.data = log2,
                   scale = TRUE,
                   cluster.idents = FALSE,
                   cluster.features = FALSE,
                   cols = c("blue","green","yellow","orange","chocolate1","red"))+
        ggtitle(group)+
        theme(axis.line=element_blank(),
                text = element_text(size=16),
                panel.grid.major = element_blank(),
                axis.text.x = element_text(size=16,
                                           angle = 90, 
                                           hjust = 0,
                                           vjust= 0.5),
                axis.title.x = element_blank(),
                axis.text.y = element_text(size=11),
                axis.title.y = element_blank(),
                plot.title = element_text(face = "plain",hjust = 0.5))+
        scale_y_discrete(position = "right")+
        coord_flip()+NoLegend()
    Progress(i, length(superfamily))
}

# for bottom 3 sub figures ==========
Enrichr_list <- list()
for(k in seq_along(csv_list)){
    csv = csv_list[k]
    Enrichr_res <- read.csv(paste0("Yang/Lung_30/GSEA/Enrichr/EVGs/",csv),stringsAsFactors = F)
    Enrichr_res = Enrichr_res[Enrichr_res$Adjusted.P.value <= adj,]
    sapply(Enrichr_res$Term, function(x) {
        Term = strsplit(x, " ")[[1]]
        paste(Term[(3-k):min(length(Term),3)],collapse = " ")
    }) %>% 
        gsub(" female","",.) %>% 
        gsub(" male","",.) -> Enrichr_res$tissue
    Enrichr_list[[k]] = Enrichr_res
}
Enrichr_df <- bind_rows(Enrichr_list)

table(Enrichr_df$tissue,Enrichr_df$cell.types) %>% 
    prop.table(2) %>% 
    as.data.frame.matrix -> prop_table
    
Enrichr_df %>%
    group_by(tissue, cell.types) %>%
    summarise(mean_Combined.Score=(mean(Combined.Score))) %>%
    spread(key = cell.types, value = mean_Combined.Score) %>% 
    column_to_rownames(var = "tissue") -> Score_table

for(i in seq_along(superfamily)){
    group = superfamily[i]
    cell.types = cell.type_list[[i]]
    cell.types = cell.types[cell.types %in% colnames(prop_table)]
    tissue = unlist(selected_tissues)
    prop_df = prop_table[tissue,cell.types]*100
    #write.csv(prop_df, file = paste0(save.path,"Percent_",names(cell.type_list)[k],"_",csv))
    Score_df = Score_table[tissue,cell.types]
    Score_df[is.na(Score_df)] =0
    #write.csv(Score_df, file = paste0(save.path,"Combined.Score_",names(cell.type_list)[k],"_",csv))
    g[[i+3]] <- DotPlot.2(Score_df,prop_df, features = cell.types, id = tissue)+NoLegend()
}

jpeg(paste0(save.path,"Dotplot.jpeg"), units="in", width=24, height=12,res=600)
wrap_plots(g, ncol = 3, nrow = 2)
dev.off()



#' @param Score_df score or expression matrix, demonstrated in color
#' @param prop_df score or expression matrix, demonstrated in dot size
#' @param features colnames of Score_df and prop_df
#' @param id rownames of Score_df and prop_df
DotPlot.2 <- function(Score_df,prop_df, features = cell.types, id = tissue, 
          log.data = log2,
          color.by = "avg.exp",
          cols = c("blue","green","yellow","orange","chocolate1","red"),
          col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
          scale = FALSE, scale.by = 'radius', split.by = NULL,split.colors = FALSE,
          scale.min = NA,scale.max = NA,exp.min = NA,exp.max = NA,n.breaks = NULL){
    prop_df %<>% rownames_to_column(var = "id") %>%
        gather("features.plot","pct.exp",-id)
    Score_df %<>% rownames_to_column(var = "id") %>%
        gather("features.plot","avg.exp",-id)
    if(identical(prop_df[,1:2],Score_df[,1:2])){
        data = cbind(prop_df, "avg.exp" = Score_df$avg.exp)
    } else stop("prop_df and Score_df have different rows or columns")
    data$id %<>% factor(levels = rev(id))
    if(is.function(log.data)) data$avg.exp = log.data(data$avg.exp+1)
    avg.exp.scaled <- sapply(
        X = unique(x = data$features.plot),
        FUN = function(x) {
            data.use <- data[data$features.plot == x, 'avg.exp']
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min, max = col.max)
            } else {
                data.use <- log(x = data.use)
            }
            return(data.use)
        }
    )
    data$avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    data$features.plot %<>% factor(levels = features)
    plot <- ggplot(data = data, mapping = aes_string(x = 'features.plot', y = 'id')) +
        geom_point(mapping = aes_string(size = 'pct.exp', fill = color.by),
                   color = "black", pch=21) +
        #scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        guides(size = guide_legend(title = 'Percent Expressed')) +
        labs(
            x = 'Features',
            y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
        ) +
        cowplot::theme_cowplot()
    
    if (split.colors) {
        plot <- plot + scale_color_identity()
    } else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    } else if (length(x = cols) == 2){
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    } else {
        plot <- plot + scale_fill_gradientn(
            colours=cols,
            n.breaks = n.breaks,
            name = ifelse(test = is.function(log.data),
                          yes = expression(atop(log[2](Combined.Score+1))),
                          no = expression(atop(Combined.Score))), 
            space = "Lab",
            na.value = "grey50",
            guide = "colourbar",
            aesthetics = "fill"
        )
    }
    
    if (!split.colors) {
        plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
    }
    plot = plot + theme(axis.line=element_blank(),
                  text = element_text(size=16),
                  panel.grid.major = element_blank(),
                  axis.text.x = element_text(size=16,
                                             angle = 90, 
                                             hjust = 1,
                                             vjust= 0.5),
                  axis.text.y = element_text(size=14),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  plot.title = element_text(face = "plain"))
    return(plot)
}
