########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","magrittr","data.table","pbapply","reshape2",
                   "tibble","tidyr","patchwork","cowplot"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

set.seed(101)
# ========== define dotplot rownames and colnames =================
cell.type_list <- list("Epithelial" = c("BC","IC","S1","S-Muc","TASC","H","p-C","C1","C-s",
                                        "Ion","NE","ME","G-Muc","G-Ser","AT1","AT2"),
                       "Structural" = c("Cr","Fb1","Fb2","Fb3","Fb4","Gli","SM1","SM2","SM3",
                                        "Pr","En-a","En-c1","En-ca","En-v","En-SM","En-l"),
                       "Immune" = c("Neu","Mon","M1","M1-2","M2","cDC","pDC","MC","B",
                                    "PC","Tcn","T-ifn","Trm","CD8-T1","T-NK","NK"))

selected_tissues <- list("GTEx tissue" = c("Lung",
                                           "Esophagus",
                                           "Skin",
                                           "Stomach",
                                           #"Salivary",
                                           "Testis",
                                           "Brain",
                                           #"Adipose",
                                           #"Vessel",
                                           "Blood",
                                           "Nerve"#,"Spleen"
                                           ),
                         "Human Gene Atlas" = c("Trachea",
                                                "Olfactory",
                                                #"HBEC",
                                                #"Tongue",
                                                #"CD33+ Myeloid",
                                                "CD14+ Mon",
                                                "BDCA4+ DC",
                                                "CD19+ B",
                                                "CD4+ T",
                                                "CD8+ T",
                                                "CD56+ NK"))
csv_list <- c("enrichR_GTEx_Tissue_Sample_Gene_Expression_Profiles_up.csv",
              "enrichR_Human_Gene_Atlas.csv")

superfamily <- c("Epithelial","Structural","Immune")
#======1.2 load  Seurat =========================
# load files
object = readRDS(file = "data/Lung_SCT_30_20210831.rds")

DefaultAssay(object) = "SCT"
meta.data = readRDS("output/20211222/meta.data_SCINA_Lung30_Azimuth_Cell_Types_2021.rds")
table(rownames(object@meta.data) == rownames(meta.data))
object@meta.data = meta.data
object$Cell_subtype %<>% factor(levels = c(unlist(cell.type_list,use.names = FALSE),"T-un","Un"))
table(object$Cell_subtype)

#object@meta.data = meta.data
object %<>% subset(subset = Doublets == "Singlet" &
                       Superfamily != "Un")
Idents(object) = "Cell_subtype"

len <- 30
#==================================
#dotplot_df <- readxl::read_excel("doc/20211006_40-gene for for dotplot revised.xlsx")
dotplot_df <- readxl::read_excel("doc/20220508_30-gene dot plot.xlsx")
dotplot_df <- dotplot_df[1:len,superfamily]
dotplot_df = dotplot_df[complete.cases(dotplot_df),]
genes = unique(unlist(df2list(dotplot_df),use.names = F))
table(genes %in% rownames(object))
# for top 3 figures ===========
g <- list()
for(i in seq_along(superfamily)){
    group = superfamily[i]
    sub_object <- subset(object, idents =  cell.type_list[[group]])
    Idents(sub_object) %<>% factor(levels = cell.type_list[[group]])
    features = dotplot_df[[group]]
    g[[i]] <- DotPlot.1(sub_object, features = rev(features),
                        log.data = log2,exp.max = 7,dot.scale = 6, 
                       scale = FALSE,
                       cluster.idents = FALSE,
                       cluster.features = FALSE,
                       cols = c("blue","green","yellow","orange","chocolate1","red"))+
        ggtitle(group)+
        theme(axis.line=element_blank(),
                text = element_text(size=16),#16
                panel.grid.major = element_blank(),
                axis.text.x = element_text(size=16,#16
                                           angle = 90, 
                                           hjust = 0,
                                           vjust= 0.5),
                axis.title.x = element_blank(),
                axis.text.y = element_text(size=14),#14
                axis.title.y = element_blank(),
                plot.title = element_text(face = "plain",hjust = 0.5))+
        scale_y_discrete(position = "right")+
        coord_flip()
    if(i < 3) g[[i]] = g[[i]] + NoLegend()
    Progress(i, length(superfamily))
}

# === generate figure =======
layout_largerFont <- c(
    area(1, 1, 14, 12),
    area(1, 13, 14, 24),
    area(1, 25, 14, 36)
)
#plot(layout_largerFont)

jpeg(paste0(path,"Dotplot_v3.jpeg"), units="in", width=19, 
     height= 8.5,res=900)
print(wrap_plots(g, nrow = 14, ncol = 36,design = layout_largerFont))
dev.off()
# for bottom 3 sub figures combine ==========
adj = 10^-4
exp.max = 300


Enrichr_list <- list()
prop_table_list <- list()
Score_table_list <- list()
save.path = path
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

for(k in 1:length(csv_list)){
    csv = csv_list[k]
    Enrichr_res <- read.csv(paste0("Yang/Lung_30/hg38/GSEA/Enrichr/",csv),stringsAsFactors = F)
    Enrichr_res = Enrichr_res[Enrichr_res$Adjusted.P.value <= adj,]
    sapply(Enrichr_res$Term, function(x) {
        Term = strsplit(x, " ")[[1]]
        paste(Term[(3-k):min(length(Term),3)],collapse = " ")
    }) %>% 
        gsub(" female","",.) %>% 
        gsub(" male","",.) -> Enrichr_res$tissue
    if(k == 1) {
        Enrichr_res$tissue %<>% tolower %>% Hmisc::capitalize()
        Enrichr_res$tissue %<>% plyr::mapvalues(from = c("Salivary gland","Blood vessel","Adipose tissue"),
                                                to   = c("Salivary", "Vessel","Adipose"))
    }
    if(k == 2) {
        Enrichr_res$tissue %<>% plyr::mapvalues(from = c("Lung",# to remove
                                                         "Testis",# to remove
                                                         "BDCA4+ DendriticCells",
                                                         "CD14+ Monocytes",
                                                         "CD19+ BCells(neg. sel.)",
                                                         "CD4+ Tcells",
                                                         "CD8+ Tcells",
                                                         "CD56+ NKCells",
                                                         "OlfactoryBulb",
                                                         "BronchialEpithelialCells"),
                                                to   = c("lung",# to remove
                                                         "testis",# to remove
                                                         "BDCA4+ DC",
                                                         "CD14+ Mon",
                                                         "CD19+ B",
                                                         "CD4+ T",
                                                         "CD8+ T",
                                                         "CD56+ NK",
                                                         "Olfactory",
                                                         "HBEC"))
    }
    Enrichr_list[[k]] = Enrichr_res

    table(Enrichr_res$tissue,Enrichr_res$cell.types) %>% 
        prop.table(2) %>% 
        as.data.frame.matrix -> prop_table_list[[k]]
    Enrichr_res %>%
        group_by(tissue, cell.types) %>%
        summarise(mean_Combined.Score=(mean(Combined.Score))) %>%
        spread(key = cell.types, value = mean_Combined.Score) %>% 
        column_to_rownames(var = "tissue") -> Score_table_list[[k]]
}

bind_rows(prop_table_list) -> prop_table
bind_rows(Score_table_list) -> Score_table


for(i in seq_along(superfamily)){
    group = superfamily[i]
    cell.types = cell.type_list[[i]]
    cell.types = cell.types[cell.types %in% colnames(prop_table)]
    empty = cell.types[!(cell.types %in% colnames(prop_table))]
    if(length(empty) >0){
        for(emp in empty) {
            prop_table %<>% cbind("empty" = 0)
            Score_table %<>% cbind("empty" = 0)
        }
        colnames(prop_table)[colnames(prop_table) %in% "empty"] = empty
        colnames(Score_table)[colnames(Score_table) %in% "empty"] = empty
    }
    tissue = unlist(selected_tissues)
    prop_table[is.na(prop_table)] =0
    prop_df = prop_table[tissue,cell.types]*100
    
    
    Score_df = Score_table[tissue,cell.types]
    Score_df[is.na(Score_df)] =0
    g[[i+3]] <- DotPlot.2(Score_df,prop_df, features = cell.types, 
                          id = tissue, scale = FALSE,log.data = NULL,
                          scale.by = "size", scale.max  = 90,
                          col.min = 0, exp.max = 300,dot.scale = 6)
    if(i < 3) g[[i+3]] = g[[i+3]] + NoLegend()
}


# === generate figure =======
layout_largerFont <- c(
    area(1, 1, 14, 12),
    area(1, 13, 14, 24),
    area(1, 25, 14, 36),
    area(15, 1, 20, 12),
    area(15, 13, 20, 24),
    area(15, 25, 20, 36)
)
#plot(layout_largerFont)

jpeg(paste0(save.path,"Dotplot_v2.jpeg"), units="in", width=19, 
     height= 14,res=900)
print(wrap_plots(g, nrow = 20, ncol = 36,design = layout_largerFont))
dev.off()

Score_table[is.na(Score_table)] = 0
prop_table[is.na(prop_table)] = 0
dotplot_res = list("Combined.Score" =Score_table[tissue,unlist(cell.type_list)],
                   "Percentage.Expressed" = prop_table[tissue,unlist(cell.type_list)]*100)
openxlsx::write.xlsx(dotplot_res, file = paste0(save.path,"Dotplot_Enricher_data.xlsx"),
                     colNames = TRUE,rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#================ Azimuth cell types ===============

adj = 2.5*10^-4
exp.max = c(300,5000)[2]
Odds.ratio = 100

Enrichr_res <- readxl::read_excel(paste0("Yang/Lung_30/hg38/GSEA/Enrichr/",
                                         "enrichR_celltypes_FC1_Azimuth_HuBMAP.xlsx"),
                                  sheet = "Azimuth_Cell_Types_2021")
#Enrichr_res = Enrichr_res[Enrichr_res$Adjusted.P.value <= adj,]
Enrichr_res = Enrichr_res[Enrichr_res$Odds.Ratio > Odds.ratio,]

Enrichr_res$tissue = gsub(" CL.*| UBER.*","",Enrichr_res$Term)
Enrichr_res$Term_CL = gsub(".* CL","CL",Enrichr_res$Term)
Enrichr_res$Term_CL = gsub(".* UBERON","UBERON",Enrichr_res$Term_CL)
Enrichr_res = Enrichr_res[order(Enrichr_res$Combined.Score,decreasing = T),]
#Enrichr_res = Enrichr_res[!duplicated(Enrichr_res$cell.types),]

table(Enrichr_res$Term,Enrichr_res$cell.types) %>% 
    prop.table(2) %>% 
    as.data.frame.matrix -> prop_table -> df
Enrichr_res %>%
    group_by(Term, cell.types) %>%
    summarise(mean_Combined.Score = mean(Combined.Score)) %>%
    spread(key = cell.types, value = mean_Combined.Score) %>% 
    column_to_rownames(var = "Term") -> Score_table
Enrichr_res %>%
    group_by(Term, cell.types) %>%
    summarise(mean_Odds.Ratio = mean(Odds.Ratio)) %>%
    spread(key = cell.types, value = mean_Odds.Ratio) %>% 
    column_to_rownames(var = "Term") -> OR_table
Enrichr_res %>%
    mutate(log.adj.p = -log10(Adjusted.P.value)) %>%
    group_by(Term, cell.types) %>%
    summarise(mean.log.adj.p = mean(log.adj.p)) %>%
    spread(key = cell.types, value = mean.log.adj.p) %>% 
    column_to_rownames(var = "Term") -> logp_table

superfamily <- c("Epithelial","Structural","Immune")
Score_table[is.na(Score_table)] = 0
logp_table[is.na(logp_table)] = 0
OR_table[is.na(OR_table)] = 0
prop_table = prop_table*100

#Score_table = sweep(Score_table, 2, colSums(Score_table),"/")
#logp_table = sweep(logp_table, 2, colSums(logp_table),"/")*100

pct.titles = c("Percent Expressed"," -log(FDR)")
score_titles = c("Combined Score", "Odds Ratio")
for(score_title in score_titles){
    save.path = paste0(path,score_title,"/")
    if(!dir.exists(save.path))dir.create(save.path, recursive = T)
    
    for(family in superfamily){
        cell.types = unlist(cell.type_list[family],use.names = F)
        table(cell.types %in% colnames(df))
        cell.types = cell.types[cell.types %in% colnames(df)]
        df1 = df[,cell.types]
        for(c in rev(cell.types)){
            df1 = df1[order(df1[,c],decreasing = T),]
        }
        
        selected_tissues = rownames(df1)[rowSums(df1) > 0]
        
        g <- DotPlot.2(Score_df = switch(score_title,
                                         "Combined Score" = Score_table,
                                         "Odds Ratio" = OR_table),
                       prop_df = logp_table,
                       features = cell.types, 
                       score_title = score_title,
                       pct.title = expression(-log[10](FDR)),
                       id = selected_tissues, scale = FALSE,log.data = NULL,
                       scale.by = "size", 
                       scale.max  = 20,
                       col.min = 0, exp.max = switch(score_title,
                                                     "Combined Score" = 5000,
                                                     "Odds Ratio" = 1000),
                       dot.scale = 6)
        
        if(family == "Epithelial"){
            jpeg(paste0(save.path,"Dotplot_Azimuth_",family,".jpeg"), units="in", 
                 width=10,height= 6,res=600)
            print(g)
            dev.off()
        }
        if(family == "Structural"){
            jpeg(paste0(save.path,"Dotplot_Azimuth_",family,".jpeg"), units="in", 
                 width=12, height= 7,res=600)
            print(g)
            dev.off()
        }
        if(family == "Immune"){
            jpeg(paste0(save.path,"Dotplot_Azimuth_",family,".jpeg"), units="in", 
                 width=10, height= 18,res=600)
            print(g)
            dev.off()
        }
    }
}



#================ 


for(k in seq_along(csv_list)){
    prop_table_list[[k]][is.na(prop_table_list[[k]])] = 0
    Colnames = cell.types[cell.types %in% colnames(prop_table_list[[k]])]
    write.csv(prop_table_list[[k]][tissue,Colnames]*100,
              file =  paste0(save.path,"Percent_",csv_list[k]))
    Score_table_list[[k]][is.na(Score_table_list[[k]])] = 0
    write.csv(Score_table_list[[k]][tissue,Colnames],
              file =  paste0(save.path,"Combined.Score_",csv_list[k]))
}

par(mfrow = c(2,1))
p1 <- unlist(prop_table_list[[1]])*100
p1 = p1[p1>25]
p2 <- unlist(prop_table_list[[2]])*100
p2 = p2[p2>25]
hist(p1,breaks = 30,xlab = "percentage",
     main = "Histogram of GTEx percentage")
hist(p2,breaks = 30,xlab = "percentage",
     main = "Histogram of HGA percentage")


par(mfrow = c(2,1))
s1 <- unlist(Score_table_list[[1]])
p1 <- unlist(prop_table_list[[1]])*100
s1 = s1[p1>=1]

s2 <- unlist(Score_table_list[[2]])
p2 <- unlist(prop_table_list[[2]])*100
s2 = s2[p1>=1]

hist(s1,breaks = 30,xlab = "Combine.score",
     main = "Histogram of GTEx Combine.score with p > 1")
hist(s2,breaks = 30,xlab = "Combine.score",
     main = "Histogram of HGA Combine.score with p > 1")


#' @param log.data = log2, 
#' @param Score_df score or expression matrix, demonstrated in color
#' @param prop_df score or expression matrix, demonstrated in dot size
#' @param features colnames of Score_df and prop_df
#' @param id rownames of Score_df and prop_df
DotPlot.2 <- function(Score_df,prop_df, features = NULL,id = NULL,
                      log.data = NULL,
                      score_title = "Combined Score", pct.title = "Percent Expressed",
                      cols = c("blue","green","yellow","orange","chocolate1","red"),
                      col.min = 0,
                      col.max = 1000,
                      dot.min = 0,
                      dot.scale = 4,
                      scale = FALSE, scale.by = 'radius', split.by = NULL,split.colors = FALSE,
                      scale.min = NA,scale.max = NA,exp.min = NA,exp.max = NA,n.breaks = NULL){
    if(is.null(features)) features = colnames(Score_df)
    if(is.null(id)) id = rownames(Score_df)
    
    prop_df %<>% rownames_to_column(var = "id") %>%
        gather("features.plot","pct.exp",-id)
    Score_df %<>% rownames_to_column(var = "id") %>%
        gather("features.plot","avg.exp",-id)
    if(identical(prop_df[,1:2],Score_df[,1:2])){
        data.plot = cbind(prop_df, "avg.exp" = Score_df$avg.exp)
    } else stop("prop_df and Score_df have different rows or columns")
    if(all(id %in% data.plot$id)){
        if(!all(data.plot$id %in% id)) data.plot = data.plot[data.plot$id %in% id,]
        data.plot$id %<>% factor(levels = rev(id))
    } else stop("some id don't exsit in rows of prop_df and Score_df")
    if(all(features %in% data.plot$features.plot)){
        if(!all(data.plot$features.plot %in% features)) data.plot = data.plot[data.plot$features.plot %in% features,]
        data.plot$features.plot %<>% factor(levels = features)
    } else stop("some id don't exsit in columns of prop_df and Score_df")
    avg.exp.scaled <- sapply(
        X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- log.data(data = data.use, min = col.min, max = col.max)
            } else if(is.function(log.data)){
                data.use <- log.data(x = data.use+1)
            }
            return(data.use)
        }
    )
    data.plot$avg.exp.scaled <- as.vector(x = avg.exp.scaled)
    color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
    color.by <- ifelse(test = is.function(log.data), yes = 'avg.exp.scaled', no = color.by)
    
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
    }
    
    if (!is.na(x = exp.min)) {
        data.plot[data.plot[,color.by] < exp.min, color.by] <- exp.min
    }
    if (!is.na(x = exp.max)) {
        data.plot[data.plot[,color.by] > exp.max, color.by] <- exp.max
    }
    scale.func <- switch(
        EXPR = scale.by,
        'size' = scale_size,
        'radius' = scale_radius,
        stop("'scale.by' must be either 'size' or 'radius'")
    )
    plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
        geom_point(mapping = aes_string(size = 'pct.exp', fill = color.by),
                   color = "black", pch=21) +
        scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        guides(size = guide_legend(title = pct.title)) +
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
                          yes = paste("log",score_title),
                          no = score_title), 
            space = "Lab",
            na.value = "grey50",
            guide = "colourbar",
            aesthetics = "fill"
        )
    }
    
    if (!split.colors) {
        plot <- plot + guides(color = guide_colorbar(title = pct.title))
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


object$cell_group = plyr::mapvalues(object$cell_types, 
                                    from = cell.type_list[["Epithelial"]],
                                    to = rep("Epithelial",length(cell.type_list[["Epithelial"]])))
object$cell_group %<>% plyr::mapvalues(from = cell.type_list[["Stromal"]],
                                       to = rep("Stromal",length(cell.type_list[["Stromal"]])))
object$cell_group %<>% plyr::mapvalues(from = cell.type_list[["Immune"]],
                                    to = rep("Immune",length(cell.type_list[["Immune"]])))
table(object$cell_group)
