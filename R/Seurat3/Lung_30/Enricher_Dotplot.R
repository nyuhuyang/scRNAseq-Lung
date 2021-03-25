########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","magrittr","data.table","pbapply",
                   "tibble","tidyr","patchwork","cowplot"), function(x) {
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
                       "Stromal"=c("F1","F2","F3","F4","Cr","Gli","Nr","SM1",
                                      "SM2","SM3","Pr","En-a","En-c","En-c1","En-v","En-l","En-sm","En-p"),
                       "Immune" = c("MC","Neu","Mon","M0","M1","M1-2","M2","M-p","DC","p-DC",
                                    "B","PC","T-cn","T-reg","T-int","T-rm","T-NK","T-ifn","T-p"))
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
csv_list <- c("enrichR_EVG_GTEx_Tissue_Sample_Gene_Expression_Profiles_up.csv",
              "enrichR_EVG_Human_Gene_Atlas.csv")
adj = c(10^-2,10^-3,10^-4,10^-5)[1]

save.path = paste0(path,"adj < ",adj,"/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

superfamily <- c("Epithelial","Stromal","Immune")
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
len <- 40

dotplot_df <- readxl::read_excel("doc/20210322_40-gene for for dotplot revised.xlsx")
dotplot_df <- dotplot_df[1:len,superfamily]
dot_features = c(dotplot_df[,superfamily[1]],
                 dotplot_df[,superfamily[2]],
                 dotplot_df[,superfamily[3]]) %>% unlist
gene_names = GTEx$V1
dot_features[!(dot_features %in% gene_names)]
all_genes = rownames(object)
dot_features[!(dot_features %in% all_genes)]
dot_features[!(dot_features %in% all_genes)]

replace_gene_names <- read.csv("Yang/GTEx/replace_gene_names.csv")
plyr::mapvalues(dot_features[!(dot_features %in% all_genes)], from = replace_gene_names$Single.cell,
                to = replace_gene_names$GTEx.v8)

newnames = plyr::mapvalues(all_genes, from = replace_gene_names$Single.cell,
                           to = replace_gene_names$GTEx.v8)
# RenameGenesSeurat  ------------------------------------------------------------------------------------
# https://github.com/satijalab/seurat/issues/2617
# Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
RenameGenesSeurat <- function(obj , newnames) { 
    print("Run this before integration. It only changes obj@assays$SCT@counts, @data and @scale.data.")
    SCT <- obj@assays$SCT
    
    if (nrow(SCT) == length(newnames)) {
        if (length(SCT@counts)) SCT@counts@Dimnames[[1]]            <- newnames
        if (length(SCT@data)) SCT@data@Dimnames[[1]]                <- newnames
        #if (length(SCT@scale.data)) SCT@scale.data@Dimnames[[1]]    <- newnames
    } else {"Unequal gene sets: nrow(SCT) != nrow(newnames)"}
    obj@assays$SCT <- SCT
    return(obj)
}
object = RenameGenesSeurat(obj = object, newnames = newnames)

# for top 3 figures ===========
g <- list()
for(i in seq_along(superfamily)){
    group = superfamily[i]
    sub_object <- subset(object, idents =  cell.type_list[[group]])
    Idents(sub_object) %<>% factor(levels = cell.type_list[[group]])
    features = features_list[[group]]
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

# for bottom 3 sub figures combine ==========
adj = 10^-4
exp.max = 300


Enrichr_list <- list()
prop_table_list <- list()
Score_table_list <- list()
save.path = path
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

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
    #cell.types = cell.types[cell.types %in% colnames(prop_table)]
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
layout <- c(
    area(1, 1, 13, 10),
    area(1, 11, 13, 18),
    area(1, 19, 13, 27),
    area(14, 1, 20, 10),
    area(14, 11, 20, 18),
    area(14, 19, 20, 27)
)
layout_largerFont <- c(
    area(1, 1, 14, 13),
    area(1, 14, 14, 24),
    area(1, 25, 14, 36),
    area(15, 1, 20, 13),
    area(15, 14, 20, 24),
    area(15, 25, 20, 36)
)
plot(layout_largerFont)

jpeg(paste0(save.path,"Dotplot_final.jpeg"), units="in", width=19, 
     height= 14,res=900)
print(wrap_plots(g, nrow = 20, ncol = 36,design = layout_largerFont))
dev.off()

for(k in seq_along(csv_list)){
    prop_table_list[[k]][is.na(prop_table_list[[k]])] = 0
    write.csv(prop_table_list[[k]]*100,
              file =  paste0(save.path,"Percent_",csv_list[k]))
    Score_table_list[[k]][is.na(Score_table_list[[k]])] = 0
    write.csv(Score_table_list[[k]],
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


#' @param Score_df score or expression matrix, demonstrated in color
#' @param prop_df score or expression matrix, demonstrated in dot size
#' @param features colnames of Score_df and prop_df
#' @param id rownames of Score_df and prop_df
DotPlot.2 <- function(Score_df,prop_df, features = cell.types, id = tissue, 
          log.data = NULL,
          cols = c("blue","green","yellow","orange","chocolate1","red"),
          col.min = 0,
          col.max = 1000,
          dot.min = 0,
          dot.scale = 4,
          scale = FALSE, scale.by = 'radius', split.by = NULL,split.colors = FALSE,
          scale.min = NA,scale.max = NA,exp.min = NA,exp.max = NA,n.breaks = NULL){
    prop_df %<>% rownames_to_column(var = "id") %>%
        gather("features.plot","pct.exp",-id)
    Score_df %<>% rownames_to_column(var = "id") %>%
        gather("features.plot","avg.exp",-id)
    if(identical(prop_df[,1:2],Score_df[,1:2])){
        data.plot = cbind(prop_df, "avg.exp" = Score_df$avg.exp)
    } else stop("prop_df and Score_df have different rows or columns")
    if(all(id %in% data.plot$id)){
        data.plot$id %<>% factor(levels = rev(id))
    } else stop("some id don't exsit in rows of prop_df and Score_df")
    if(all(features %in% data.plot$features.plot)){
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
                          yes = expression(atop("Mean enrichR's\n Combined.Score",log[2](Combined.Score+1))),
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
