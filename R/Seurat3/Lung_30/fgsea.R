# https://github.com/ctlab/fgsea/releases/tag/v1.15.1
devtools::install_github("ctlab/fgsea")
library(dplyr)
library(magrittr)
library(stringr)
library(ggsci)
library(fgsea)
library(openxlsx)
library(progress)


path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

# prepare gmt file
repare_gmt_file = FALSE
if(repare_gmt_file) {
    read.path = "../seurat_resources/msigdb/"
    folder_list <- c("@1-Development","@2-Physiology","@3-Disease","@4-Cancer")
    for(i in seq_along(folder_list)){
        pathway_name_list <- list.files(paste0(read.path,folder_list[i])) %>% 
            gsub("\\.txt","",.) %>% gsub("^[0-9][0-9]_","",.)
        txt_list <- list.files(paste0(read.path,folder_list[i]), full.names = T)
        pathway <- lapply(txt_list, function(x) {
            read.table(x, sep = "\t",stringsAsFactors = F) %>% 
                pull %>% gsub(" ///.*","",.) %>% gsub(",.*","",.)})
        names(pathway) = pathway_name_list
        fgsea::writeGmtPathways(pathway, paste0(read.path,sub(".*-","",folder_list[i]),
                                                ".gmt"))
    }
}

# read pathway
hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v7.1.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))

Biocarta <- gmtPathways("../seurat_resources/msigdb/c2.cp.biocarta.v7.1.symbols.gmt")
kegg <- gmtPathways("../seurat_resources/msigdb/c2.cp.kegg.v7.1.symbols.gmt")
pid <- gmtPathways("../seurat_resources/msigdb/c2.cp.pid.v7.1.symbols.gmt")
reactome <- gmtPathways("../seurat_resources/msigdb/c2.cp.reactome.v7.1.symbols.gmt")
tft <- gmtPathways("../seurat_resources/msigdb/c3.tft.v7.1.symbols.gmt")
Development <- gmtPathways("../seurat_resources/msigdb/Development.gmt")
Physiology <- gmtPathways("../seurat_resources/msigdb/Physiology.gmt")
Disease <- gmtPathways("../seurat_resources/msigdb/Disease.gmt")
Cancer <- gmtPathways("../seurat_resources/msigdb/Cancer.gmt")

msigdb_list <- list("hallmark" = hallmark,
                    "Biocarta" = Biocarta,
                    "kegg" = kegg,
                    "pid" = pid,
                    "reactome" = reactome,
                    "transcription_factor_targets" = tft,
                    "Development" = Development,
                    "Physiology" = Physiology,
                    "Disease" = Disease,
                    "Cancer" = Cancer)

set.seed(101)
read.path = "Yang/Lung_30/DE_analysis/"
#===========================
# read data
# C_Cell_types
project = "C_Cell_types/"
save.path = paste0("Yang/Lung_30/GSEA/","fgsea/",project)
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

csv_list <- list.files(paste0(read.path, project), full.names = T)
deg_list <- lapply(csv_list, function(x) read.csv(x, row.names = 1,
                                                  stringsAsFactors=F))
deg_list %<>% lapply(function(res) res[order(res["p_val_adj"]),])
(clusters <- sapply(deg_list, function(res) unique(res$cluster)))
names(deg_list) = clusters
res =  bind_rows(deg_list)
(clusters = unique(as.character(res$cluster)))

fgseaRes <- list()
pb <- progress_bar$new(total = length(clusters))
for(i in 1:length(clusters)){
    message(clusters[i])
    geneRank = res[res$cluster == clusters[i],]
    #geneRank = geneRank[geneRank$avg_logFC >= 0.5 & geneRank$p_val_adj < 0.05, ]
    geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")]  %>% tibble::deframe()
    
    tmp <- list()
    for(m in 1:length(msigdb_list)){
        #print(paste(m,names(msigdb_list)[m]))
        tmp[[m]] <- fgseaMultilevel(pathways=msigdb_list[[m]], 
                                    stats=geneRank, scoreType = "pos",eps = 0,
                                    BPPARAM=SerialParam())
    }

    # record and remove empty element in tmp
    emp <- c()
    for(k in seq_along(msigdb_list)) {
        if(!is.null(tmp[[k]]) ) {
            tmp[[k]][, library := names(msigdb_list[k])]
        } else emp = c(emp, k)
    }
    if(!is.null(emp)) tmp[emp] = NULL
    
    tmp = bind_rows(tmp)
    tmp = tmp[tmp$padj < 0.05,]
    if(nrow(tmp) > 0 ) {
        tmp$cluster = clusters[i]
    } else tmp =NULL
    fgseaRes[[i]] = tmp
    pb$tick()
}

df_fgseaRes =  bind_rows(fgseaRes)
df_fgseaRes = df_fgseaRes[df_fgseaRes$padj < 0.05, ]

fgseaRes_list <- split(df_fgseaRes, f = df_fgseaRes$library)
for(i in seq_along(fgseaRes_list)){
    fgseaRes_list[[i]]$leadingEdge =NULL
    write.csv(fgseaRes_list[[i]], file = paste0(save.path,"fgsea_celltypes_",
                                                   names(fgseaRes_list)[i],".csv"),
              row.names = FALSE)
}
#========= A_Sample_types ==================
set.seed(101)
read.path = "Yang/Lung_30/DE_analysis/"
project = "A_Sample_types/"
save.path = paste0("Yang/Lung_30/GSEA/","fgsea/",project)
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

csv_list <- list.files(paste0(read.path, project), full.names = T)
deg_list <- list()
for(i in seq_along(csv_list)){
    if(file.size(csv_list[i]) > 10) deg_list[[i]] = read.csv(csv_list[i],
                                                             row.names = 1,
                                                             stringsAsFactors=F)
    Progress(i, length(csv_list))
}
Names <- gsub(".*Lung_30_A_","",csv_list) %>% 
    gsub("^.*_celltypes=","",.) %>%
    gsub("\\.csv","",.)
tail(Names)
names(deg_list) = Names
deg_list = deg_list[lapply(deg_list,length)>0]
Names = names(deg_list)
deg_list %<>% lapply(function(res) res[order(res["p_val_adj"]),])
(clusters <- lapply(deg_list, function(res) unique(res$cluster)))

# select colored cells

wb <- loadWorkbook("doc/20200903_Comparison groups of cells - Yang modified.xlsx")
df_annotation <- readxl::read_excel("doc/20200903_Comparison groups of cells - Yang modified.xlsx",
                                    sheet = "A - Sample types")
# filter by color
df_color <- data.frame()
for(i in 1:length(wb$styleObjects)){
    x = wb$styleObjects[[i]]
    if(is.null(x$style$fill)) next
    print(i)
    if(x$style$fill$fillFg == "FFFFFF00" & x$sheet == "A - Sample types"){
        message(i)
        df_color = rbind(df_color, data.frame(ri =x$rows, ci = x$cols, col = "FFFFFF00"))
    }
}
df_color$ri %<>% -1
select_groups = df_annotation[unique(df_color$ri), c(1:2)] %>% as.data.frame()
colnames(select_groups)[1] = "k"
select_groups$k %<>% as.integer()
select_groups$`CELL GROUPS` %<>% gsub(" \\(.*", "", .)
select_groups = plyr::arrange(select_groups, k) # r order numeric values
group_names = sapply(1:23, function(x) paste0(select_groups$k[x],"_", select_groups$`CELL GROUPS`[x]))
group_names %<>% gsub("\\+","\\\\+",.)
print(selected_groups <- grep(paste(group_names,collapse = "|"), Names, value = T) %>%
          sort)

fgseadRes <- list()
pb <- progress_bar$new(total = length(selected_groups))
for(i in 1:length(selected_groups)){
    res <- deg_list[[selected_groups[i]]]
    tem_list <- list()
    for(m in 1:2){
        geneRank = res[res$cluster == clusters[[selected_groups[i]]][m],]
        geneRank = geneRank[log10(geneRank$p_val_adj) <= -20, ]
        geneRank = geneRank[geneRank$avg_logFC > 0, ]
        geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")] %>% tibble::deframe()
        if(length(geneRank) < 10) next
        tmp <- lapply(msigdb_list, function(x) fgseaMultilevel(pathways=x, 
                                                               stats=geneRank, 
                                                               scoreType = "pos",
                                                               eps = 0,
                                                               BPPARAM=BiocParallel::SerialParam())) # for linux
        # record and remove empty element in tmp
        emp <- c()
        for(k in seq_along(msigdb_list)) {
            if(!is.null(tmp[[k]]) ) {
                tmp[[k]][, library := names(msigdb_list[k])]
            } else emp = c(emp, k)
        }
        if(!is.null(emp)) tmp[emp] = NULL
        
        tmp = bind_rows(tmp)
        tmp = tmp[tmp$padj < 0.05,]
        if(nrow(tmp) > 0 ) {
            tmp$cluster = clusters[[selected_groups[i]]][m]
        } else tmp =NULL
        tem_list[[m]] = tmp
    }
    fgseadRes[[i]] = bind_rows(tem_list)
    fgseadRes[[i]]$groups = selected_groups[i]
    pb$tick()
}

for(i in 1:length(fgseadRes)){
    fgseadRes[[i]]$leadingEdge = NULL
    Progress(i, length(fgseadRes))
    }
    
df_fgseadRes =  bind_rows(fgseadRes)
df_fgseadRes = df_fgseadRes[df_fgseadRes$padj < 0.05, ]

fgseadRes_list <- split(df_fgseadRes, f = df_fgseadRes$library)
for(i in seq_along(fgseadRes_list)){
    write.csv(fgseadRes_list[[i]], file = paste0("Yang/Lung_30/GSEA/fgsea/A_Sample_types/fgsea_",
                                                   names(fgseadRes_list)[i],".csv"),
              row.names = FALSE)
    Progress(i, length(fgseadRes_list))
}

#========= B_Cell_groups ==================
project = "B_Cell_groups/"
save.path = paste0("Yang/Lung_30/GSEA/fgsea/",project)
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

csv_list <- list.files(paste0(read.path, project), full.names = T)
deg_list <- list()
for(i in seq_along(csv_list)){
    if(file.size(csv_list[i]) > 10) deg_list[[i]] = read.csv(csv_list[i],
                                                             row.names = 1,
                                                             stringsAsFactors=F)
    Progress(i, length(csv_list))
}
Names <- gsub(".*Lung_30_B_","",csv_list) %>% 
    gsub("\\.csv","",.)
tail(Names)
names(deg_list) = Names
deg_list = deg_list[lapply(deg_list,length)>0]
Names = names(deg_list)
deg_list %<>% lapply(function(res) res[order(res["p_val_adj"]),])
(clusters <- lapply(deg_list, function(res) unique(res$cluster)))

# select colored cells

wb <- loadWorkbook("doc/20200903_Comparison groups of cells - Yang modified.xlsx")
df_annotation <- readxl::read_excel("doc/20200903_Comparison groups of cells - Yang modified.xlsx",
                                    sheet = "B - Cell groups")
# filter by color
df_color <- data.frame()
for(i in 1:length(wb$styleObjects)){
    x = wb$styleObjects[[i]]
    if(is.null(x$style$fill)) next
    print(i)
    if(x$style$fill$fillFg == "FFFFFF00" & x$sheet == "B - Cell groups"){
        message(i)
        df_color = rbind(df_color, data.frame(ri =x$rows, ci = x$cols, col = "FFFFFF00"))
    }
}
df_color$ri %<>% -1
select_groups <- data.frame()
for(m in 1:nrow(df_color)) {
    k =  df_annotation[df_color$ri[m], 1] %>% pull
    group = df_annotation[df_color$ri[m], df_color$ci[m]] %>% pull
    select_groups = rbind(select_groups ,data.frame("k" = k,
                                                    "CELL.GROUPS" = group,
                                                    stringsAsFactors = F))
    
}

colnames(select_groups)[1] = "k"
select_groups$k %<>% as.integer()
select_groups$`CELL.GROUPS` %<>% gsub(" \\(.*", "", .)
select_groups %<>% plyr::arrange(k) # r order numeric values
group_names = apply(select_groups,1, function(x) paste0(x["k"],"_", x["CELL.GROUPS"]))
group_idx <- gsub("_.*","",group_names) %>% as.integer() %>% unique() %>% sort
select_files <- data.frame(k = group_idx, csv_files = Names[group_idx],stringsAsFactors = F)
df_select_groups <- full_join(select_groups, select_files, by = "k")

fgseadRes <- list()
pb <- progress_bar$new(total = nrow(df_select_groups))
for(i in 1:nrow(df_select_groups)){
    res <- deg_list[[df_select_groups$csv_files[i]]]
    geneRank = res[res$cluster == df_select_groups$CELL.GROUPS[i],]
    #geneRank = geneRank[log10(geneRank$p_val_adj) <= -20, ]
    geneRank = geneRank[geneRank$avg_logFC > 0, ]
    geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")]  %>% tibble::deframe()
    print(paste("gene number = ",length(geneRank)))
    if(length(geneRank) < 2) next
    tmp <- lapply(msigdb_list, function(x) fgseaMultilevel(pathways=x, 
                                                           stats=geneRank, 
                                                           scoreType = "pos",
                                                           eps = 0))
    # record and remove empty element in tmp
    emp <- c()
    for(k in seq_along(msigdb_list)) {
        if(!is.null(tmp[[k]]) ) {
            tmp[[k]][, library := names(msigdb_list[k])]
        } else emp = c(emp, k)
    }
    if(!is.null(emp)) tmp[emp] = NULL
    tmp = bind_rows(tmp)
    tmp = tmp[tmp$padj < 0.05,]
    if(nrow(tmp) > 0 ) {
        tmp$cluster = df_select_groups$CELL.GROUPS[i]
        tmp$groups = sub("^([0-9][0-9])_","",df_select_groups$csv_files[i])
    } else tmp =NULL
    
    fgseadRes[[i]] = tmp
    pb$tick()
}

df_fgseadRes =  bind_rows(fgseadRes)
df_fgseadRes$leadingEdge = NULL
df_fgseadRes = df_fgseadRes[df_fgseadRes$padj < 0.05, ]

fgseadRes_list <- split(df_fgseadRes, f = df_fgseadRes$library)
for(i in seq_along(fgseadRes_list)){
    write.csv(fgseadRes_list[[i]], file = paste0("Yang/Lung_30/GSEA/fgsea/B_Cell_groups/fgsea_",
                                                   names(fgseadRes_list)[i],".csv"),
              row.names = FALSE)
    Progress(i, length(fgseadRes_list))
}