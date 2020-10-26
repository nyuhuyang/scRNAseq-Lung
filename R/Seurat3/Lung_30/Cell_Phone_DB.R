conda create --name python3.6.7 python=3.6.7
conda activate python3.6.7
conda install -c conda-forge r-base r-ggplot2 r-pheatmap r-dplyr r-tidyr
pip install cellphonedb

#conda activate cellphonedb
cd Yang/Lung_30/Cell_Phone_DB/
    
    for con in proximal distal terminal COPD:
    mkdir $con
cd $con
done
cellphonedb method statistical_analysis ../Lung_30-${con}_meta.data_cellphone.txt ../Lung_30-${con}_counts.txt
cellphonedb plot dot_plot
cellphonedb plot heatmap_plot ../Lung_30-${con}_meta.data_cellphone.txt

########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
####################################
invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#===== seperate by cell types ======================================================================
read.path = "Yang/Lung_30/Cell_Phone_DB/"
# prepare label matrix
sections = c("proximal","distal","terminal","COPD")
labels = c("AT1","AT2","AT2-1","AT2-p","B","BC","BC-p","BC-S",
           "C1","C2","C3","Cr","DC","En-A","En-C","En-C1","En-L",
           "En-p","En-SM","En-V","F1","F2","F3","F4","Gli","H","IC-S",
           "IC1","IC2","Ion","M-p","M0","M1","M1-2","M2","MC","MEC",
           "Mon","NEC","Neu","Nr","p-C","P-DC","PC","Pr","RBC","S",
           "S-d","SM1","SM2","SM3","SMG-Muc","SMG-Ser","T-cn","T-ifn",
           "T-int","T-NK","T-p","T-reg","T-rm","T-un","T7","Un")
N = length(labels)
labels_mt <- matrix(paste0(rep(labels, each = N),"|",
                                 rep(labels, times = N)),
                    nrow = N,ncol = N,
                    dimnames = list(labels,labels))
# declare temp results list
txt_list <- vector("list", length = N*4)
names(txt_list) =  paste0(labels,".",rep(sections,each = N))

# load data
for(i in seq_along(sections)){
    section = sections[i]
    means <- readxl::read_excel(paste0(read.path,section,"/out/means.xlsx"))
    pvalues <- readxl::read_excel(paste0(read.path,section,"/out/pvalues.xlsx"))
    for(m in seq_along(labels)) {
        label.section = paste0(labels[m],".",section)
        label.available = labels_mt[,labels[m]] %in% colnames(means) %>% labels_mt[,labels[m]][.]
        n = length(label.available)
        if(n > 0) {
            txt_list[[label.section]] = cbind(means[,1:11], means[,label.available])
            colnames(txt_list[[label.section]])[12:(12+n-1)] %<>% paste0(".",section,".means")
            txt_list[[label.section]] %<>% cbind(pvalues[,label.available])
            colnames(txt_list[[label.section]])[(12+n):(12+2*n-1)] %<>% paste0(".",section,".pvalues")
        }
    }
    svMisc::progress(i/length(sections)*100)
}

label.section.remove = sapply(txt_list,is.null)
txt_list = txt_list[!label.section.remove]

# declare results list
results_list <- vector("list", length = N)
names(results_list) =  labels

# sepcify columns to joyn
join_by_columns = c("id_cp_interaction","interacting_pair","partner_a","partner_b","gene_a","gene_b",
                "secreted","receptor_a","receptor_b","annotation_strategy","is_integrin")
# sepcify custom order to reorder the columns. 
# Total should have 63 labels * 63 labels* 4 conditions * 2 values = 31752 combinantions
category.all = paste0(rep(as.vector(labels_mt),each = 8),".",
                      rep(rep(sections,each = 2),times = 63*63),".",
                      rep(c("means","pvalues"),times = 63*63*4))

for(m in seq_along(labels)) {
    label = labels[m]
    label.section.avaible = grep(label,names(txt_list),value = T)
    tmp = txt_list[label.section.avaible] %>% 
        purrr::reduce(dplyr::full_join, by = join_by_columns)
    
    # remove all 0 in mean and all 1 in pavlue
    tmp_mean = data.matrix(tmp[,grep("mean",colnames(tmp),value = T)],rownames.force = TRUE)
    tmp_mean[tmp_mean == 0] <- ""
    
    tmp_pvalue = data.matrix(tmp[,grep("pvalue",colnames(tmp),value = T)],rownames.force = TRUE)
    tmp_pvalue[tmp_pvalue == 1] <- ""
    
    #Reorder rows using custom order
    tmp = cbind(tmp[,join_by_columns],tmp_mean)
    tmp %<>% cbind(tmp_pvalue)
    
    category.available = category.all[category.all %in% colnames(tmp)]
    tmp = tmp[,c(join_by_columns, category.available)]

    results_list[[label]] = tmp
    svMisc::progress(m/length(labels)*100)
}

for(i in seq_along(results_list)){
    openxlsx::write.xlsx(results_list[i], 
                         file = paste0(path,"Cell_Phone_DB_",names(results_list)[i],".xlsx"),
                         colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
    svMisc::progress(i/length(results_list)*100)
}

for(i in seq_along(results_list)){
    write.table(results_list[[i]],file = paste0(path,names(results_list)[i],".txt"),
                sep = "\t",row.names = FALSE)
    svMisc::progress(i/length(results_list)*100)
}


#===== seperate by regions ======================================================================
read.path = "Yang/Lung_30/Cell_Phone_DB/"
# prepare label matrix
sections = c("proximal","distal","terminal","COPD")
sections_inital = toupper(substr(sections,1,1))
# sepcify columns to joyn
join_by_columns = c("id_cp_interaction","interacting_pair","partner_a","partner_b","gene_a","gene_b",
                    "secreted","receptor_a","receptor_b","annotation_strategy","is_integrin","receptor_a|b")

# load data
mean_pvalues_list <- vector("list", length = 4)
names(mean_pvalues_list) = toupper(substr(sections,1,1))
for(i in seq_along(sections)){
    section = sections[i]
    means <- readxl::read_excel(paste0(read.path,section,"/out/means.xlsx"))
    pvalues <- readxl::read_excel(paste0(read.path,section,"/out/pvalues.xlsx"))
    means %<>% gather(key = "receptor_a|b", value = "Mean",
                       -colnames(means)[1:11])
    pvalues %<>%  gather(key = "receptor_a|b", value = "p value",
                       -colnames(means)[1:11])
    
    mean_pvalues = dplyr::full_join(means,pvalues, by = join_by_columns)
    mean_pvalues = mean_pvalues[mean_pvalues$`p value` <= 0.05,]
    mean_pvalues = mean_pvalues[order(mean_pvalues$Mean,decreasing = TRUE),]
    colnames(mean_pvalues)[13:14] %<>% paste(toupper(substr(section,1,1)))
    mean_pvalues_list[[i]] = mean_pvalues
    
    svMisc::progress(i/length(sections)*100)
}

# load data
full_mean_pvalues_list <- vector("list", length = 4)
names(full_mean_pvalues_list) = toupper(substr(sections,1,1))

for(i in seq_along(sections)){

    # rearrange mean_pvalues_list
    reOrder_sections_inital = c(toupper(substr(sections[i],1,1)),
                                sections_inital[-i])
    mean_pvalues_list = mean_pvalues_list[reOrder_sections_inital]
    full_mean_pvalues_list[[i]] = mean_pvalues_list %>% purrr::reduce(dplyr::left_join, by = join_by_columns)
    svMisc::progress(i/length(sections)*100)
}


for(i in seq_along(full_mean_pvalues_list)){
    openxlsx::write.xlsx(full_mean_pvalues_list[i], 
                         file = paste0(path,"Cell_Phone_DB_",names(full_mean_pvalues_list)[i],".xlsx"),
                         colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
    svMisc::progress(i/length(full_mean_pvalues_list)*100)
}

openxlsx::write.xlsx(full_mean_pvalues_list, 
                     file = paste0(path,"Cell_Phone_DB_by_regions.xlsx"),
                     colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
