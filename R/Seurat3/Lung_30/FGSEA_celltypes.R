library(dplyr)
library(magrittr)
library(stringr)
library(ggsci)
# https://github.com/ctlab/fgsea/releases/tag/v1.15.1
install.packages("~/Downloads/fgsea-1.15.0.tar.gz", repos=NULL, type="source")


library(fgsea)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

set.seed(101)
read.path = "Yang/Lung_30/DE_analysis/"

# read pathway
hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v7.1.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))

Biocarta <- gmtPathways("../seurat_resources/msigdb/c2.cp.biocarta.v7.1.symbols.gmt")
kegg <- gmtPathways("../seurat_resources/msigdb/c2.cp.kegg.v7.1.symbols.gmt")
pid <- gmtPathways("../seurat_resources/msigdb/c2.cp.pid.v7.1.symbols.gmt")
reactome <- gmtPathways("../seurat_resources/msigdb/c2.cp.reactome.v7.1.symbols.gmt")
tft <- gmtPathways("../seurat_resources/msigdb/c3.tft.v7.1.symbols.gmt")

msigdb_list <- list("hallmark" = hallmark,
                    "Biocarta" = Biocarta,
                    "kegg" = kegg,
                    "pid" = pid,
                    "reactome" = reactome,
                    "transcription_factor_targets" = tft)

#===========================
# read data
project = "C_Cell_types/"
save.path = paste0(path,project)
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

csv_list <- list.files(paste0(read.path, project), full.names = T)
deg_list <- lapply(csv_list, function(x) read.csv(x, row.names = 1,
                                                  stringsAsFactors=F))
deg_list %<>% lapply(function(res) res[order(res["p_val_adj"]),])
(clusters <- sapply(deg_list, function(res) unique(res$cluster)))
names(deg_list) = clusters


groups <- list("EPITHELIAL"="BC+BC-S+BC-p+IC1+IC2+IC-S+S+S-d+p-C+H+C1+C2+C3+SMG-Muc+SMG-Ser+Ion+NEC+MEC",
               "STROMAL" ="F1+F2+F3+F4+SM1+SM2+SM3+Pr+Cr+Gli",
               "ENDOTHELIAL" = "En-A+En-C+En-C1+En-V+En-L+En-p",
               "T_CELLS" = "T-cn+T-rm+T7+T-NK+T-ifn+T-int+T-p+T-reg+T-un",
               "IMMUNE" = "Mon+M0+M1+M2+DC+M1-2+M-p+MC+Neu+P-DC+B+PC+T-cn+T-rm+T7+T-NK+T-ifn+T-reg+T-p+T-int+T-un")
for(m in seq_along(msigdb_list)){
    save.path.sub = paste0(save.path,names(msigdb_list)[m],"/")
    if(!dir.exists(save.path.sub))dir.create(save.path.sub, recursive = T)
    
    for(i in seq_along(groups)) {
        name = names(groups[i])
        group = groups[[i]]
        group = str_split(group,pattern = "\\+") %>% unlist()
        res = bind_rows(deg_list[group])
        fgseaRes = FgseaDotPlot(stats=res, pathways = msigdb_list[[m]],
                                padj = 0.25,pval = 0.05,
                                decreasing = F,
                                Rowv = T,Colv = T,
                                size = " -log10(pval)", fill = "NES",
                                pathway.name = names(msigdb_list)[m],rotate.x.text = T,
                                title = paste(names(msigdb_list)[m],"enriched in",name),
                                font.xtickslab=12, font.main=12, font.ytickslab = 10,
                                font.legend = list(size = 12),font.label = list(size = 12),
                                do.return = T,save.path = save.path.sub, do.print = T,
                                width = length(group)/5+6,height = length(msigdb_list[[m]])/50+5,hjust = 0.5)
        write.csv(fgseaRes, file = paste0(save.path.sub,"fgsea_celltypes_",
                                          names(msigdb_list)[m],"_",name,"_FDR0.25_pval0.05.csv"))
    }
    Progress(m, length(msigdb_list))
}


# combine =====

for(m in 6:length(msigdb_list)){
    message(names(msigdb_list)[m])
    res =  bind_rows(deg_list)

    (clusters = unique(as.character(res$cluster)))
    fgseaRes <- list()
    for(i in 1:length(clusters)){
        geneRank = res[res$cluster == clusters[i],]
        geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")]  %>% tibble::deframe()
        fgseaRes[[i]] <- fgseaMultilevel(pathways=msigdb_list[[m]], stats=geneRank,nproc = 1)
        fgseaRes[[i]] = as.data.frame(fgseaRes[[i]])
        fgseaRes[[i]] = fgseaRes[[i]][,c("pathway","pval","padj","NES")]
        fgseaRes[[i]] = fgseaRes[[i]][fgseaRes[[i]]$pval < 0.05,]
        fgseaRes[[i]] = fgseaRes[[i]][fgseaRes[[i]]$padj < 0.25,]
        if(nrow(fgseaRes[[i]]) > 0 ) {
            fgseaRes[[i]]$cluster = clusters[i]
        } else fgseaRes[[i]] =NULL
        Progress(i, length(clusters))
    }
    df_fgseaRes <- bind_rows(fgseaRes)
    df_fgseaRes = df_fgseaRes[!is.na(df_fgseaRes[, "pathway"]),]
    df_fgseaRes[," -log10(pval)"] = -log10(df_fgseaRes$pval)
    df_fgseaRes[," -log10(padj)"] = -log10(df_fgseaRes$padj)
    write.csv(df_fgseaRes, file = paste0(save.path,"fgsea_celltypes_",
                                      names(msigdb_list)[m],"_FDR0.25_pval0.05.csv"),
              row.names = FALSE)
}

#===========================
project = "A_Sample_types/"
save.path = paste0(path,project)
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

csv_list <- list.files(paste0(read.path, project), full.names = T)
deg_list <- list()
for(i in seq_along(csv_list)){
    if(file.size(csv_list[i]) > 10) deg_list[[i]] = read.csv(csv_list[i],
                                                             row.names = 1,
                                                             stringsAsFactors=F)
    Progress(i, length(csv_list))
}

for(k in seq_along(deg_list)){
    stats = deg_list[[k]]
    clusters = unique(as.character(stats$cluster))
    fgseaRes <- list()
    for(i in seq_along(clusters)){
        geneRank = stats[stats$cluster == clusters[i],]
        geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")]  %>% tibble::deframe()
        fgseaRes[[i]] <- fgseaMultilevel(pathways=hallmark, stats=geneRank)
        fgseaRes[[i]] = as.data.frame(fgseaRes[[i]])
        fgseaRes[[i]] = fgseaRes[[i]][,c("pathway","pval","padj","NES")]
        if(clusters[i] == order.yaxis.by[1]) {
            order.yaxis = fgseaRes[[i]][order(fgseaRes[[i]][,order.yaxis.by[2]],
                                              decreasing = decreasing), "pathway"]
        }
        if(!is.null(pval)) fgseaRes[[i]] = fgseaRes[[i]][fgseaRes[[i]]$pval < pval,]
        if(!is.null(padj)) fgseaRes[[i]] = fgseaRes[[i]][fgseaRes[[i]]$padj < padj,]
        if(nrow(fgseaRes[[i]]) > 0 ) {
            fgseaRes[[i]]$cluster = clusters[i]
        } else fgseaRes[[i]] =NULL
        Progress(i, length(clusters))
    }
    df_fgseaRes <- data.table::rbindlist(fgseaRes) %>% as.data.frame()
    df_fgseaRes = df_fgseaRes[!is.na(df_fgseaRes[, "pathway"]),]
    df_fgseaRes[," -log10(pval)"] = -log10(df_fgseaRes$pval)
    df_fgseaRes[," -log10(padj)"] = -log10(df_fgseaRes$padj)
}
    

xlsx_color_theme(path, check_filetype = TRUE)
geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")]  %>% tibble::deframe()
fgseaRes[[i]] <- fgseaMultilevel(pathways=pathways, stats=geneRank)