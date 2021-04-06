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
object <- readRDS(file = "data/Lung_SCT_30_20200710.rds") 

DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
Idents(object) = "cell_types"
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

