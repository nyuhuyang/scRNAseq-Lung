# conda activate r4.0
library(dplyr)
library(magrittr)
library(stringr)
library(ggsci)
library(enrichR)
library(openxlsx)
library(progress)
library(pbapply)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

set.seed(101)
read.path = "Yang/Lung_30/DE_analysis/"

# select pathway
wb <- loadWorkbook("Yang/Lung_30/GSEA/Enrichr/Enrichr libraries.xlsx")
df <- readxl::read_excel("Yang/Lung_30/GSEA/Enrichr/Enrichr libraries.xlsx")

# filter by color
df_color <- data.frame()
for(i in 1:length(wb$styleObjects)){
    x = wb$styleObjects[[i]]
    if(is.null(x$style$fill)) next
    if(x$style$fill$fillFg == "FFFFFF00"){
        df_color = rbind(df_color, data.frame(ri =x$rows, ci = x$cols, col = "FFFFFF00"))
    }
}
dbs = df[unique(df_color$ri)-1, "Gene-set Library"] %>% pull %>% sort
dbs = dbs[-which(dbs %in% "dbGaP")]

#===========================
project = "A_Sample_types/"
save.path = paste0(path,"Enrichr/",project)
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

enrichedRes <- list()
pb <- progress_bar$new(total = length(selected_groups))
for(i in 1:length(selected_groups)){
    res <- deg_list[[selected_groups[i]]]
    tem_list <- list()
    for(m in 1:2){
        geneRank = res[res$cluster == clusters[[selected_groups[i]]][m],]
        #geneRank = geneRank[log10(geneRank$p_val_adj) <= -20, ]
        geneRank = geneRank[geneRank$avg_logFC > 0, ]
        geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")]  %>% tibble::deframe()
        print(paste("gene number = ",length(geneRank)))
        if(length(geneRank) < 10) next
        invisible(capture.output(tmp <- enrichr(names(geneRank), dbs)))
        # record and remove empty element in tmp
        emp <- c()
        for(k in seq_along(tmp)) {
            if(nrow(tmp[[k]]) > 0 ) {
                tmp[[k]][,"library"] = names(tmp[k])
            } else emp = c(emp, k)
        }
        if(!is.null(emp)) tmp[emp] = NULL
        
        tmp = bind_rows(tmp)
        tmp = tmp[tmp$Adjusted.P.value < 0.05,]
        if(nrow(tmp) > 0 ) {
            tmp$cluster = clusters[[selected_groups[i]]][m]
            tmp$groups = selected_groups[i]
        } else tmp =NULL
        tem_list[[m]] = tmp
    }
        enrichedRes[[i]] = bind_rows(tem_list)
        #write.csv(enrichedRes[[i]], file = paste0(save.path,"enrichR_",
        #                             selected_groups[i],".csv"),
        #          row.names = FALSE)
    pb$tick()
}

df_enrichedRes =  bind_rows(enrichedRes)
df_enrichedRes = df_enrichedRes[df_enrichedRes$Adjusted.P.value < 0.05, ]

enrichedRes_list <- split(df_enrichedRes, f = df_enrichedRes$library)
for(i in seq_along(enrichedRes_list)){
    #enrichedRes_list[[i]] = enrichedRes_list[[i]][,-grep("X.html.|X.",colnames(enrichedRes_list[[i]]))]
    write.csv(enrichedRes_list[[i]], file = paste0("Yang/Lung_30/GSEA/Enrichr/A_Sample_types/enrichR_",
                                                   names(enrichedRes_list)[i],".csv"),
              row.names = FALSE)
    Progress(i, length(enrichedRes_list))
}
#===========================
project = "B_Cell_groups/"
save.path = paste0("Yang/Lung_30/GSEA/Enrichr/",project)
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

enrichedRes <- list()
pb <- progress_bar$new(total = nrow(df_select_groups))
for(i in 1:nrow(df_select_groups)){
    res <- deg_list[[df_select_groups$csv_files[i]]]

    geneRank = res[res$cluster == df_select_groups$CELL.GROUPS[i],]
    geneRank = geneRank[log10(geneRank$p_val_adj) <= -20, ]
    geneRank = geneRank[geneRank$avg_logFC > 0, ]
    geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")]  %>% tibble::deframe()
    print(paste("gene number = ",length(geneRank)))
    if(length(geneRank) < 2) next
    tmp <- enrichr(names(geneRank), dbs)
    # record and remove empty element in tmp
    emp <- c()
    for(k in seq_along(tmp)) {
        if(nrow(tmp[[k]]) > 0 ) {
            tmp[[k]][,"library"] = names(tmp[k])
        } else emp = c(emp, k)
    }
    if(!is.null(emp)) tmp[emp] = NULL
    
    tmp = bind_rows(tmp)
    tmp = tmp[tmp$Adjusted.P.value < 0.05,]
    if(nrow(tmp) > 0 ) {
        tmp$cluster = df_select_groups$CELL.GROUPS[i]
        tmp$groups = sub("^([0-9][0-9])_","",df_select_groups$csv_files[i])
    } else tmp =NULL

    enrichedRes[[i]] = tmp
    #write.csv(enrichedRes[[i]], file = paste0(save.path,"enrichR_",
    #                                          df_select_groups$csv_files[i],".csv"),
    #          row.names = FALSE)
    pb$tick()
}

df_enrichedRes =  bind_rows(enrichedRes)
df_enrichedRes = df_enrichedRes[df_enrichedRes$Adjusted.P.value < 0.05, ]

enrichedRes_list <- split(df_enrichedRes, f = df_enrichedRes$library)
for(i in seq_along(enrichedRes_list)){
    write.csv(enrichedRes_list[[i]], file = paste0("Yang/Lung_30/GSEA/Enrichr/B_Cell_groups/enrichR_",
                                                   names(enrichedRes_list)[i],".csv"),
              row.names = FALSE)
    Progress(i, length(enrichedRes_list))
}

#===========================
# read data
# C_Cell_types
project = "C_Cell_types/"
save.path = paste0(path,"Enrichr/",project)
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

csv_list <- list.files(paste0(read.path, project), full.names = T)
deg_list <- lapply(csv_list, function(x) read.csv(x, row.names = 1,
                                                  stringsAsFactors=F))
deg_list %<>% lapply(function(res) res[order(res["p_val_adj"]),])
(clusters <- sapply(deg_list, function(res) unique(res$cluster)))
names(deg_list) = clusters
res =  bind_rows(deg_list)
(clusters = unique(as.character(res$cluster)))

enrichedRes <- list()
for(i in 1:length(clusters)){
    geneRank = res[res$cluster == clusters[i],]
    geneRank = geneRank[geneRank$avg_logFC >= 0.5 & geneRank$p_val_adj < 0.05, ]
    geneRank = geneRank[order(geneRank["avg_logFC"]),c("gene","avg_logFC")]  %>% tibble::deframe()
    tmp <- enrichr(names(geneRank), dbs) #dbs[-which(dbs %in% "dbGaP")])
    # record and remove empty element in tmp
    emp <- c()
    for(k in seq_along(tmp)) {
        if(nrow(tmp[[k]]) > 0 ) {
            tmp[[k]][,"library"] = names(tmp[k])
        } else emp = c(emp, k)
    }
    if(!is.null(emp)) tmp[emp] = NULL
    
    tmp = bind_rows(tmp)
    tmp = tmp[tmp$Adjusted.P.value < 0.05,]
    if(nrow(tmp) > 0 ) {
        tmp$cluster = clusters[i]
    } else tmp =NULL
    enrichedRes[[i]] = tmp
    Progress(i, length(clusters))
}
df_enrichedRes =  bind_rows(enrichedRes)
df_enrichedRes = df_enrichedRes[df_enrichedRes$Adjusted.P.value < 0.05, ]
write.xlsx(df_enrichedRes, asTable = F,
           file = paste0(save.path,"enrichR_celltypes_FDR0.05.xlsx"),
           borders = "surrounding")

enrichedRes_list <- split(df_enrichedRes, f = df_enrichedRes$library)
for(i in seq_along(enrichedRes_list)){
    write.csv(enrichedRes_list[[i]], file = paste0(save.path,"enrichR_celltypes_",
                                                   names(enrichedRes_list)[i],".csv"),
              row.names = FALSE)
    Progress(i, length(enrichedRes_list))
    
}
write.xlsx(enrichedRes_list, asTable = F,
           file = paste0(save.path,"enrichR_celltypes_FDR0.25_pval0.05.xlsx"),
           borders = "surrounding")

##################
projects = c("C_Cell_types/","A_Sample_types/adj p < 10(-20)",
            #"B_Cell_groups/adj p < 10(-20)",
            "B_Cell_groups/adj p < 0.05")
for(k in seq_along(projects)){
    save.path = paste0("Yang/Lung_30/GSEA/Enrichr/",projects[k])
    csv_list <- list.files(save.path, full.names = T)
    deg_list <- lapply(csv_list, function(x) read.csv(x, #row.names = 1,
                                                      stringsAsFactors=F))
    deg_list <- lapply(deg_list, function(x) {
        x$Overlap = paste0(" ",x$Overlap)
        return(x)}
    )
    for(i in seq_along(deg_list)){
        write.csv(deg_list[[i]], file = csv_list[i],row.names = FALSE)
        Progress(i, length(deg_list))
    }
}


#====cell type-specific EVG gene ==============
project = "EVGs/"
save.path = paste0("Yang/Lung_30/GSEA/Enrichr/",project)
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

superfamily <- c("Epithelial","Structural","Immune")
res_list <- pblapply(superfamily, function(sheet) {
    readxl::read_excel("Yang/Lung_30/DE_analysis/F_EVGs_allCells/Lung_30-EVGs-full.xlsx",
                       sheet = sheet)
    })
cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix("", n-nrow(x), ncol(x))))) 
}

res_df <- do.call(cbind.fill, res_list)

df <- readxl::read_excel("doc/Chord diagram cell order - updated abbreviations 12-14-20.xlsx",col_names = T)
(cell.types <- sort(df$cell.types))
deg_list <- pblapply(cell.types, function(cell){
    tmp = res_df[,grep(paste0("^",cell,"_"),colnames(res_df),value = T)] %>% as.data.frame()
    tmp = tmp[!is.na(tmp[,paste0(cell,"_gene")]),]
    tmp = tmp[tmp[,paste0(cell,"_gene")] != "",]
    tmp[,2:6] %<>% apply(2,as.numeric)
    tmp
})
names(deg_list) = cell.types
deg_list = deg_list[lapply(deg_list,length)>0]

enrichedRes <- list()
for(i in 1:length(cell.types)){
    cell = cell.types[i]
    geneRank =  deg_list[[cell.types[i]]]
    colnames(geneRank) %<>% gsub(paste0(cell,"_"),"",.)
    geneRank = geneRank[order(geneRank["FC"]),c("gene","FC")]  %>% tibble::deframe()
    tmp <- enrichr(names(geneRank), dbs) #dbs[-which(dbs %in% "dbGaP")])
    # record and remove empty element in tmp
    emp <- c()
    for(k in seq_along(tmp)) {
        if(nrow(tmp[[k]]) > 0 ) {
            tmp[[k]][,"library"] = names(tmp[k])
        } else emp = c(emp, k)
    }
    if(!is.null(emp)) tmp[emp] = NULL
    
    tmp = bind_rows(tmp)
    tmp = tmp[tmp$Adjusted.P.value < 0.05,]
    if(nrow(tmp) > 0 ) {
        tmp$cell.types = cell.types[i]
    } else tmp =NULL
    enrichedRes[[i]] = tmp
    Progress(i, length(cell.types))
}
df_enrichedRes =  bind_rows(enrichedRes)
df_enrichedRes = df_enrichedRes[df_enrichedRes$Adjusted.P.value < 0.05, ]
write.xlsx(df_enrichedRes, asTable = F,
           file = paste0(save.path,"enrichR_EVG_all.xlsx"),
           borders = "surrounding")

enrichedRes_list <- split(df_enrichedRes, f = df_enrichedRes$library)
for(i in seq_along(enrichedRes_list)){
    write.csv(enrichedRes_list[[i]], file = paste0(save.path,"enrichR_EVG_",
                                                   names(enrichedRes_list)[i],".csv"),
              row.names = FALSE)
    Progress(i, length(enrichedRes_list))
    
}
write.xlsx(enrichedRes_list, asTable = F,
           file = paste0(save.path,"enrichR_EVG_celltypes.xlsx"),
           borders = "surrounding")
#how many samples (and how many different tissues) were used in this HGA data set? 
library(data.table)
HGA <- fread('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Human_Gene_Atlas',
             fill=TRUE,comment.char=.)
# For GTEx, was the v8 release of GTEx included in Enrichr? This will give an idea how many samples were included overall.
GTEx <- fread('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GTEx_Tissue_Sample_Gene_Expression_Profiles_up',
              fill=TRUE)
GTEx <- fread("/Users/yah2014/Downloads/GTEx_Tissue_Sample_Gene_Expression_Profiles_up.txt",fill=TRUE,
              sep="\t")
