invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "tidyverse","ggpubr","rstatix","kableExtra"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
save.path = "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/figures/"
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

set.seed(101)
"
Analysis 1
3 groups: D-young; D-old; COPD

Sub-analysis for all T:
        calculate “max normal level (MNL)”: mean + 2SD for IFNG expression in all T cells of D-young
proportions of cells (among all T cells) with expression > MNL in each group
output: graph; table (cell %) – per group and per sample

Sub-analysis for T-rm:
        calculate “max normal level (MNL)”: mean + 2SD for IFNG expression in T-rm cells of D-young
proportions of cells (among T-rm cells) with expression > MNL in each group
output: graph; table – per group and per sample

Sub-analysis for T-NK:
        calculate “max normal level (MNL)”: mean + 2SD for IFNG expression in T-NK cells of D-young
proportions of cells (among T-NK cells) with expression > MNL in each group
output: graph; table – per group and per sample

Analysis 2 (you did during our zoom meeting)
2 groups: D-all; COPD

Sub-analysis for all T:
        calculate “max normal level (MNL)”: mean + 2SD for IFNG expression in all T cells of D-all
proportions of cells (among all T cells) with expression > MNL in each group
output: graph; table – per group and per sample

Sub-analysis for T-rm:
        calculate “max normal level (MNL)”: mean + 2SD for IFNG expression in T-rm cells of D-all
proportions of cells (among T-rm cells) with expression > MNL in each group
output: graph; table – per group and per sample

Sub-analysis for T-NK:
        calculate “max normal level (MNL)”: mean + 2SD for IFNG expression in T-NK cells of D-all
proportions of cells (among T-NK cells) with expression > MNL in each group
output: graph; table – per group and per sample


4 - Correlation between IFNG and CD8A, IFNG and CD8B 

Single-cell data:
- D+COPD samples; all T cells only 
- D+COPD samples; T-rm only 
- D samples; all T cells only
- D-samples; T-rm only
- COPD samples; all T cells only
- COPD samples; T-rm only

GTEx data:
- all 577 GTEx samples
- all 577 GTEx samples (dots for younger colored green; dots for older colored red)
- separate for younger samples (20-29 + 30-39 combined) 
- separate for older samples (60-69 + 70-79 combined)

x-axis - IFNG; y -axis - CD8A or CD8B
Please include statistics in each plot: Spearman Rho and p value.

"
# 1 - IFNG expression in T cells (all T cell combined) of D vs COPD
T_cells = c("T-cn","T-ifn","T-int","T-NK","T-p","T-reg","T-rm")
age_df = data.frame("orig.ident" = c("UNC-54-D", "UNC-57-D", "UNC-66-D", "UNC-70-D",
                                     "CU-12-D","CU-12-D-R",
                                     "UNC-44-D", "UNC-48-D", "UNC-55-D", "UNC-67-D", "UNC-69-D", "UNC-71-D", "VU-27-D"),
                    "age" = c(rep("old",4),rep("middle",2),rep("young",7)))

object = readRDS(file = "data/Lung_SCT_30_20200710.rds")
DefaultAssay(object) = "SCT"
object %<>% subset(subset = cell_types %in% T_cells & 
                           conditions %in% c("COPD", "distal"))
table(object$orig.ident, object$cell_types) %>% kable() %>% kable_styling()
# Analysis 1
# 3 groups: D-young; D-old; COPD
object$age <- plyr::mapvalues(object$orig.ident, from = age_df$orig.ident, to =age_df$age)
object@meta.data[!object$age %in% c("old","middle","young"),"age"] = ""
object$conditions %<>% gsub("distal","D-",.)
object$conditions %<>% paste0(object$age)
table(object$conditions)
data = FetchData(object, vars = c("IFNG","conditions","cell_types","orig.ident"))
data %<>% filter(conditions %in% c("D-young","D-old","COPD")) 
data$conditions %<>% factor(levels = c("D-young","D-old","COPD")) 

Analysis_1 <- function(data, MNL = NULL,calculate.MNL.by = "D-young", filter.by){
        data %<>% filter(cell_types %in% filter.by)
        if(!is.numeric(MNL)) {
                if(!is.null(calculate.MNL.by)) {
                print(MNL <- data %>% filter(conditions %in% calculate.MNL.by) %>% 
                              summarise(MNL = mean(IFNG)+1*sd(IFNG)))
                MNL = pull(MNL)
                } else MNL = 0
        } 
        per_group <- data %>% group_by(conditions) %>%
                summarise(group = unique(conditions),
                          total_num = length(IFNG),
                          total_mean = mean(IFNG),
                          MNL = MNL,
                          high_proportion = sum(IFNG > MNL)/length(IFNG),
                          high_num = sum(IFNG > MNL),
                          high_mean = ifelse(any(IFNG > MNL),mean(IFNG),0), 
                          low_proportion = sum(IFNG <= MNL)/length(IFNG),
                          low_num = sum(IFNG <= MNL),
                          low_mean = ifelse(any(IFNG <= MNL),mean(IFNG),0))
        colnames(per_group)[1] = "sample"
        per_sample <- data %>% group_by(orig.ident) %>%
                summarise(total_num = length(IFNG),
                          total_mean = mean(IFNG),
                          MNL = MNL,
                          high_proportion = sum(IFNG > MNL)/length(IFNG),
                          high_num = sum(IFNG > MNL),
                          high_mean = ifelse(any(IFNG > MNL),mean(IFNG),0), 
                          low_proportion = sum(IFNG <= MNL)/length(IFNG),
                          low_num = sum(IFNG <= MNL),
                          low_mean = ifelse(any(IFNG <= MNL),mean(IFNG),0))
        colnames(per_sample)[1] = "sample"
        data1 = data[!duplicated(data$orig.ident),]
        per_sample$group = plyr::mapvalues(per_sample$sample,
                                           from = data1$orig.ident, to = as.character(data1$conditions))
        per_sample  %<>% arrange(desc(group))
        per_sample = per_sample[,c("sample","group",colnames(per_sample)[-c(1,ncol(per_sample))])]
        
        df = rbind(per_group, per_sample)
        return(df)
}

results1 <- pbapply::pblapply(T_cells,function(x) Analysis_1(data,filter.by = x, 
                                                             MNL = 0.7613702))
names(results1) = T_cells
results1[["all_T"]] = Analysis_1(data = data,calculate.MNL.by = "D-young", filter.by = T_cells)
results1 = results1[c("all_T",T_cells)]
openxlsx::write.xlsx(results1, file =  "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/analysis_1.xlsx",
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

data$conditions %<>% factor(levels=c("D-young", "D-old","COPD"))
Analysis_plot <- function(data, MNL = NULL, calculate.MNL.by = "D-young", label.y = NULL,
                          file.name, fun = ggviolin,x = "group", y = "high_mean",
                          formula = high_mean ~ group...){
        if(!is.numeric(MNL)) {
                if(!is.null(calculate.MNL.by)) {
                        print(MNL <- data %>% filter(conditions %in% calculate.MNL.by) %>% 
                                      summarise(MNL = mean(IFNG)+1*sd(IFNG)))
                        MNL = pull(MNL)
                } else MNL = 0
        }
        lvls = levels(data$conditions)
        per_sample <- data %>% group_by(orig.ident) %>%
                summarise(high_proportion = sum(IFNG > MNL)/length(IFNG),
                          high_num = sum(IFNG > MNL),
                          high_mean = ifelse(any(IFNG > MNL),mean(IFNG),0), 
                          low_proportion = sum(IFNG <= MNL)/length(IFNG),
                          low_num = sum(IFNG <= MNL),
                          low_mean = ifelse(any(IFNG <= MNL),mean(IFNG),0))
        colnames(per_sample)[1] = "sample"
        data1 = data[!duplicated(data$orig.ident),]
        per_sample$group = plyr::mapvalues(per_sample$sample,
                                           from = data1$orig.ident, to = as.character(data1$conditions))
        per_sample$group %<>% factor(levels = lvls)
        if(is.null(label.y)) {
                label.y = max(per_sample[,y])
                label.y = seq(from = label.y+0.2, to = label.y+0.4,by = 0.1)
        }
        my_comparisons <- list( lvls[1:2], lvls[2:3], lvls[c(1,3)] )
        jpeg(filename = file.name, units="in", width=3.5, height=3.5,res=600)
        print(fun(per_sample, x = x, y = y,
                 color = x,
                 palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                 add = c("boxplot"),error.plot = "errorbar",bxp.errorbar = TRUE,
                 bxp.errorbar.width = 0.2)+ 
                stat_compare_means(comparisons = my_comparisons, label.y = label.y)+ # Add pairwise comparisons p-value
                stat_compare_means(label.y = max(label.y)+0.1))     # Add global p-value
        dev.off()
        res = wilcox_test(formula = formula, data = per_sample)
        res %<>% cbind(kruskal_test(formula = formula, data = per_sample))
        res$MNL = MNL
        res[2:3,10:16] = ""
        return(res)
}

for(test in c("high_mean","high_proportion")){
        formula = switch (test,
                          "high_mean" = high_mean ~ group,
                          "high_proportion" = high_proportion ~ group)
        Kruskal <- list()
        Kruskal[["T_cells"]] <- Analysis_plot(data, calculate.MNL.by = "D-young",#label.y = c(0.5,0.6,0.7),
                                              fun = ggviolin,formula = formula,
                                              y = test,
                                              file.name = paste0(save.path,"scRNA-seq/",test,"/violint_D-young_D-old_COPD_T_cells.jpeg"))
        for(cell in T_cells){
                Kruskal[[cell]] = Analysis_plot(data[data$cell_types %in% cell,],MNL = 0.7613702, #label.y = c(0.7,0.6,0.5),
                                                fun = ggviolin,formula = formula,
                                                y = test,
                                                file.name = paste0(save.path,"scRNA-seq/",test,"/violin_D-young_D-old_COPD_",cell,".jpeg"))
        }
        Kruskal[["T_cells"]] <- Analysis_plot(data, calculate.MNL.by = "D-young",#label.y = c(0.5,0.6,0.7),
                                              fun = ggboxplot,formula = formula,
                                              y = test,
                                              file.name = paste0(save.path,"scRNA-seq/",test,"/boxplot_D-young_D-old_COPD_T_cells.jpeg"))
        for(cell in T_cells){
                Kruskal[[cell]] = Analysis_plot(data[data$cell_types %in% cell,],MNL = 0.7613702, #label.y = c(0.7,0.6,0.5),
                                                fun = ggboxplot,formula = formula,
                                                y = test,
                                                file.name = paste0(save.path,"scRNA-seq/",test,"/boxplot_D-young_D-old_COPD_",cell,".jpeg"))
        }
        openxlsx::write.xlsx(Kruskal, file = paste0(save.path,"scRNA-seq/",test,"/",test," p-values.xlsx"))

}

Idents(object) = "conditions"
object %<>% subset(idents = c("COPD","D-old","D-young"))
object@meta.data$conditions %<>% factor(levels = c("D-young","D-old","COPD"))




data = FetchData(object, vars = c("IFNG","CD8A","CD8B", "conditions","cell_types","orig.ident"))
data %<>% filter(conditions %in% c("COPD","D-old","D-young"))
cor_res <- Hmisc::rcorr(x = as.matrix(data[,c("IFNG","CD8A","CD8B")]),type="spearman")
df = data.frame("correlation" = cor_res$r[1,c("CD8A","CD8B")],
            "pvalue" = cor_res$P[1,c("CD8A","CD8B")],
            "markers" = c("CD8A","CD8B"))
data.table::fwrite(data, file = paste0(save.path,"scRNA-seq/Expression data used for scRNA-seq IFNG Scatter plot.csv"))
data.table::fwrite(df, file = paste0(save.path,"scRNA-seq/correlation data used for scRNA-seq IFNG Scatter plot.csv"))


jpeg(paste0(save.path,"scRNA-seq/4 -Single-cell- Correlation (scatter plots) IFNG and CD8A.jpeg"),
     units="in", width=5, height=5,res=600)
ggscatter(data, x = "IFNG", y = "CD8B", color = "conditions",
          palette = c("#228B22","#E7B800",'#601b3f','#3D9970','#FF4136','#FF851B'))+
        annotate("text", x = 3, y = 0.5, 
                 label = paste("correlation = ",round(df["CD8A","correlation"],digits = 3),
                               ", p value = ", round(df["CD8A","pvalue"],digits = 3)),
                 size = 5)
dev.off()

jpeg(paste0(save.path,"scRNA-seq/4 -scRNA- Correlation (scatter plots) IFNG and CD8A.jpeg"),
     units="in", width=5.5, height=4.5,res=600)
FeatureScatter(object,feature1 = "IFNG",feature2 = "CD8A",
               cols = c("#228B22","#E7B800",'#601b3f'))+ggtitle("")+
        annotate("text", x = 3, y = 2.5, 
                 label = paste("correlation = ",round(df["CD8A","correlation"],digits = 3),
                               ", p value = ", round(df["CD8A","pvalue"],digits = 3)),
                 size = 5)
dev.off()
jpeg(paste0(save.path,"scRNA-seq/4 -scRNA- Correlation (scatter plots) IFNG and CD8B.jpeg"),
     units="in", width=5.5, height=4.5,res=600)
FeatureScatter(object,feature1 = "IFNG",feature2 = "CD8B",
               cols = c("#228B22","#E7B800",'#601b3f'))+ggtitle("")+
        annotate("text", x = 3, y = 2.5, 
                 label = paste("correlation = ",round(df["CD8B","correlation"],digits = 3),
                               ", p value = ", round(df["CD8B","pvalue"],digits = 3)),
                 size = 5)
dev.off()

# Analysis 2 (you did during our zoom meeting)
# 2 groups: D-all; COPD
data = FetchData(object, vars = c("IFNG","conditions","cell_types","orig.ident"))
data$conditions %<>% factor(levels = c("D-young","D-old","D-middle","COPD"))
data$cell_types %<>% factor(levels = T_cells)
data = data[order(data$conditions,data$cell_types),]
data.table::fwrite(data, file = paste0(save.path,"scRNA-seq/Expression data used for scRNA-seq IFNG vilplot.csv"))

write.csv(data, file = "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/IFNG_exp_all-T.csv")
write.csv(data[data$cell_types %in% "T-rm",], file = "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/IFNG_exp_T-rm.csv")

data$conditions %<>% gsub("D-.*","D-all",.)
data.table::fwrite(data, file = paste0(save.path,"scRNA-seq/Expression data used for scRNA-seq IFNG Scatter plot.csv"))

results2 <- pbapply::pblapply(T_cells,function(x) Analysis_1(data,filter.by = x, 
                                                             MNL = 0.8980386))
names(results2) = T_cells
results2[["all_T"]] = Analysis_1(data = data, calculate.MNL.by = "D-all",filter.by = T_cells)
results2 = results2[c("all_T",T_cells)]
openxlsx::write.xlsx(results2, file =  "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/analysis_2.xlsx",
                     colNames = TRUE,row.names = T,borders = "surrounding",colWidths = c(NA, "auto", "auto"))


Analysis_plot(data, calculate.MNL.by = "D-all",
              filter.by = T_cells,ylim = c(0, 8),label.y = 6,
              file.name = "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/ggviolin plots/D-all_COPD_T_cells~.jpeg")
Analysis_plot(data[data$cell_types %in% "T-rm",],calculate.MNL.by = "D-all", 
              filter.by = T_cells,ylim = c(0, 8),label.y = 6,
              file.name = "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/ggviolin plots/D-all_COPD_T-rm~.jpeg")
Analysis_plot(data[data$cell_types %in% "T-NK",],calculate.MNL.by = "D-all", 
              filter.by = T_cells,ylim = c(0, 7),label.y = 5,
              file.name = "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/ggviolin plots/D-all_COPD_T-NK~.jpeg")
Idents(object) = "cell_types"


# 2 - IFNG expression in T cells (all T cell combined) of D-old vs D-young

data = FetchData(object, vars = c("IFNG","orig.ident","cell_types","conditions"))
data$age <- plyr::mapvalues(data$orig.ident, from = age_df$orig.ident, to =age_df$age)
data %<>% filter(cell_types %in% T_cells) %>%
        filter(age %in% c("young", "old"))
stat.test <- data %>%
        t_test(IFNG ~ age) %>%
        adjust_pvalue()
stat.test

jpeg(paste0(save.path,"scRNA-seq/2 - IFNG expression in T cells (all T cell combined) of D-old vs D-young.jpeg"),
     units="in", width=5, height=3.5,res=600)

ggviolin(data, x = "age", y = "IFNG",
         color = "age",ylim = c(0, 5.5),
         palette = c("#00AFBB", "#E7B800"))+#,add = c("boxplot"))+
        stat_pvalue_manual(stat.test, 
                           y.position = 5.0,
                           label = "T-test, p value = {p}"
        )
dev.off()

#3 - IFNG expression in GTEx samples: old (60-69 + 70-79 combined) vs young (20-29 + 30-39 combined) 
age_df = data.frame("Age" =c(20,  30,  40,  50,  60,  70),
                    "age" = c("young", "young","middle","middle","old","old"))
load("data/Lung_GTEx_20210226.Rda")
data = FetchData(object, vars = c("IFNG","Sex","Age"))
data$age <- plyr::mapvalues(data$Age, from = age_df$Age, to =age_df$age)
data %<>% filter(age %in% c("young", "old"))
data$age %<>% factor(levels = c("young","old"))
write.csv(data, file = paste0(save.path,"GTEx/Expression data used for GTEx IFNG violin plot.csv"))

(stat.test1 <- data %>%
        t_test(IFNG ~ age) %>%
        adjust_pvalue())

jpeg(paste0(save.path,"GTEx/3 - IFNG expression in GTEx samples: young (20-39) vs old (60-79)_all.jpeg"),
     units="in", width=5, height=3.5,res=600)
ggviolin(data, x = "age", y = "IFNG",
         color = "age",ylim = c(0, 9.5),
         palette = c("#228B22","#E7B800"),add = c("boxplot"))+
        stat_pvalue_manual(stat.test1, 
                           y.position = 9,
                           label = "T-test, p value = {p}"
        )
dev.off()

jpeg(paste0(save.path,"GTEx/3 - IFNG expression in GTEx samples: young (20-39) vs old (60-79)_males.jpeg"),
     units="in", width=5, height=3.5,res=600)
(stat.test2 <- data %>% filter(Sex == "male") %>%
        t_test(IFNG ~ age) %>%
        adjust_pvalue())
ggviolin(data %>% filter(Sex == "male"), x = "age", y = "IFNG",
         color = "age",ylim = c(0, 9.5),
         palette = c("#228B22","#E7B800"),add = c("boxplot"))+
        stat_pvalue_manual(stat.test2, 
                           y.position = 9,
                           label = "T-test, p value = {p}"
        )
dev.off()
jpeg(paste0(save.path,"GTEx/3 - IFNG expression in GTEx samples: young (20-39) vs old (60-79)_female.jpeg"),
     units="in", width=5, height=3.5,res=600)
(stat.test3 <- data %>% filter(Sex == "female") %>%
                t_test(IFNG ~ age) %>%
                adjust_pvalue())
ggviolin(data %>% filter(Sex == "female"), x = "age", y = "IFNG",
         color = "age",ylim = c(0, 9.5),
         palette = c("#228B22","#E7B800"),add = c("boxplot"))+
        stat_pvalue_manual(stat.test3, 
                           y.position = 9,
                           label = "T-test, p value = {p}"
        )
dev.off()

stat.test <- bind_rows(list(stat.test1,stat.test2,stat.test3))
stat.test$catalog = c("all","male","female")
write.csv(stat.test, file = paste0(save.path,"GTEx/pvalue data used for GTEx IFNG violin plot.csv"))

# 4 - Correlation between IFNG and CD8A, IFNG and CD8B 
age_df = data.frame("Age" =c(20,  30,  40,  50,  60,  70),
                    "age" = c("young", "young","middle","middle","old","old"))
load("data/Lung_GTEx_20210226.Rda")
data = FetchData(object, vars = c("IFNG","CD8A","CD8B","Sex","Age"))
data$age <- plyr::mapvalues(data$Age, from = age_df$Age, to =age_df$age)
data$age %<>% factor(levels = c("young","middle","old"))
write.csv(data, paste0(save.path,"GTEx/Expression data used for GTEx IFNG scatter plot.csv"))

GTEx_list <- pbapply::pblapply(list(data,
                                   data %>% filter(age == "young"),
                                   data %>% filter(age == "old") ), FUN = function(x){
        cor_res <- Hmisc::rcorr(x = as.matrix(x[,c("IFNG","CD8A","CD8B")]),type="spearman")
        df = data.frame("correlation" = cor_res$r[1,c("CD8A","CD8B")],
                        "pvalue" = cor_res$P[1,c("CD8A","CD8B")],
                        "markers" = c("CD8A","CD8B"))
        return(df)
})
GTEx_res = bind_rows(GTEx_list)
GTEx_res$catalog = rep(c("all","younger","older"),each =2)
write.csv(GTEx_res, file = paste0(save.path,"GTEx/Correlation data used for GTEx IFNG scatter plot.csv"))

#GTEx_res$` -log10(pvalue)` = -log10(GTEx_res$pvalue)
#GTEx_res[GTEx_res$` -log10(pvalue)` > 15, " -log10(pvalue)"] = 15
my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")#定义色板

jpeg(paste0(path,"4 -GTEx- Correlation (scatter plots) IFNG and CD8A.jpeg"),
     units="in", width=5, height=5,res=600)
ggscatter(data[data$age %in% c("old","young"),], x = "IFNG", y = "CD8A", color = "age",
          palette = c("#228B22","#E7B800"))+
        annotate("text", x = 3, y = 1, 
                 label = paste("correlation = ",round(GTEx_res["CD8A...1","correlation"],digits = 3),
                               "p value = ", GTEx_res["CD8A...1","pvalue"]),
                 size = 5)
dev.off()
jpeg(paste0(path,"4 -GTEx- Correlation (scatter plots) IFNG and CD8B.jpeg"),
     units="in", width=5, height=5,res=600)
ggscatter(data[data$age %in% c("old","young"),], x = "IFNG", y = "CD8B", color = "age",
          palette = c("#228B22","#E7B800"))+
        annotate("text", x = 3, y = 0.5, 
                 label = paste("correlation = ",round(GTEx_res["CD8B...2","correlation"],digits = 3),
                               ", p value = ", GTEx_res["CD8B...2","pvalue"]),
                 size = 5)
dev.off()

jpeg(paste0(path,"4 -GTEx- Correlation between IFNG and CD8A, IFNG and CD8B.jpeg"),
     units="in", width=5, height=5,res=600)
ggdotchart(GTEx_res, x = "markers", y ="catalog",color = "correlation",
           size = " -log10(pvalue)",rotate = T,
           ggtheme = theme_bw())
dev.off()

"""
Violin plots:

IFNG expression in immune cell types (order as in our dot-plot you generated for single-cell markers) 
IFNGR1 expression in surface airway epithelial cell types (order as listed in my email below for GSEA)
IFNGR2 expression in surface airway epithelial cell types (order as listed in my email below for GSEA)

4 variants (v1-v4) for each plot:

v1 - for all samples not divided into groups other than cell types (highlight cell types with significant enrichment; adj p<0.05 vs other cell types) 
v2 - grouped per P and D side-by-side for each cell type listed above (highlight cell types with significant difference; adj p<0.05 P vs D) 
v3 - grouped per D and COPD for each cell type listed above (highlight cell types with significant difference; adj p<0.05 D vs COPD) 
v4 - grouped per D-young and D-old for each cell type listed above (highlight cell types with significant difference; adj p<0.05 D-young vs D-old) 

Excel data: statistics for each of the 4 variants above: 

average UMI, Pct1 % positive cells, total number of cells of a given cell type 
log2FC and adj p value for comparisons:

v1 - compared to other cell types in the family shown in the plot 
v2 - between P and D
v3 - between COPD and D
v4 - between D-old and D-young
"""

object = readRDS(file = "data/Lung_30_20200710.rds")
SAE_cells = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","d-S",
           "H","p-C","C1","C2","C3","Ion","NE")
Immune_cells = c("MC","Neu","Mon","M0","M1","M1-2","M2","M-p","DC","p-DC",
           "B","PC","T-cn","T-reg","T-int","T-rm","T-NK","T-ifn","T-p")
age_df = data.frame("orig.ident" = c("UNC-54-D", "UNC-57-D", "UNC-66-D", "UNC-70-D",
                                     "UNC-44-D", "UNC-48-D", "UNC-55-D", "UNC-67-D", "UNC-69-D", "UNC-71-D", "VU-27-D"),
                    "age" = c(rep("older",4),rep("younger",7)))

Idents(object) = "cell_types"
Immune <- subset(object, idents = Immune_cells)
Immune$cell_types %<>% factor(levels = Immune_cells)
Immune %<>% ScaleData("IFNG")
SAE <- subset(object, idents = SAE_cells)
SAE$cell_types %<>% factor(levels = SAE_cells)
SAE %<>% ScaleData(features = c("IFNGR1","IFNGR2"))
# ============= Immune ==================
jpeg(paste0(path,"IFNG_Immune_v1.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(Immune, features = "IFNG",group.by = "cell_types",
        assay = "RNA",pt.size = 0.01)  + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##################  Immune VlnPlot ################## 
Idents(Immune) = "conditions"
P_D <- subset(Immune,idents = c("proximal","distal"))
P_D$conditions %<>% factor(levels = c("proximal","distal"))
D_COPD <- subset(Immune,idents = c("distal","COPD"))
D_COPD$conditions %<>% factor(levels = c("distal","COPD"))
D <- subset(Immune,idents = "distal")
D$age <- plyr::mapvalues(D$orig.ident, from = age_df$orig.ident, to =age_df$age)
Idents(D) = "age"
D %<>% subset(idents = c("younger","older"))
D$age %<>% factor(levels = c("younger","older"))

jpeg(paste0(path,"IFNG_Immune_v2.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(P_D, group.by = "cell_types",features = "IFNG",
        split.by = "conditions",assay = "RNA",pt.size = 0.01)  + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

jpeg(paste0(path,"IFNG_Immune_v3.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(D_COPD, group.by = "cell_types",features = "IFNG",
        split.by = "conditions",assay = "RNA",pt.size = 0.01)  + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

jpeg(paste0(path,"IFNG_Immune_v4.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(D, group.by = "cell_types",features = "IFNG",split.by = "age",
        assay = "RNA",pt.size = 0.01)  + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

########### Immune average expression for all ##############
AverageExpression_pct_v1 <- function(object, group.by = "cell_types", 
                                  features = "IFNG", assays = "SCT"){
        UMI_df <- FetchData(object = object, vars = c(features,group.by),)
        UMI_df$IFNG %<>% expm1
        UMI <- group_by(UMI_df,cell_types) %>% 
                summarise_at(vars(features),list(average_UMI = mean,
                                                 p_positive_cells = function(x) sum(x >0)/length(x),
                                                 n_positive_cells = function(x) sum(x >0)))
        UMI$average_UMI %<>% log1p()
        rownames(UMI) = UMI$cell_types
        return(UMI)
}

Idents(Immune) = "conditions"
P <- subset(Immune,idents = c("proximal"))
D <- subset(Immune,idents = c("distal"))
COPD <- subset(Immune,idents = c("COPD"))
D$age <- plyr::mapvalues(D$orig.ident, from = age_df$orig.ident, to =age_df$age)
Idents(D) = "age"
younger<- subset(D,idents = "younger")
older<- subset(D,idents = "older")
res_list <- list()
UMI_list1 <- pbapply::pblapply(list(Immune,P,D,younger,older,COPD), 
                              FUN = AverageExpression_pct_v1,features = "IFNG")
res_list[["IFNG immune"]] = bind_cols(UMI_list1)
########### Immune average expression for per sample ##############
orig_idents = c("UNC-48-P","UNC-55-P","UNC-66-P","UNC-69-P","UNC-71-P",
                "UNC-48-D","UNC-55-D","UNC-66-D","UNC-69-D","UNC-71-D",
                "UNC-44-D","UNC-54-D","UNC-57-D","UNC-67-D","UNC-70-D",
                "CU-12-D","CU-12-D-R","VU-27-D","UNC-48-T","UNC-55-T",
                "UNC-66-T","UNC-69-T","UNC-71-T","UNC-44-T","CU-12-T",
                "UNC-51-D","UNC-52-D","UNC-61-D","VU19-D","VU-37-D")
T_cells <- c("T-cn","T-reg","T-int","T-rm","T-NK","T-ifn","T-p")
AverageExpression_pct_v2 <- function(orig.ident = "UNC-48-P",
                                     object, 
                                     group.by = "cell_types",
                                     group.by.list = list("all immune cells" = Immune_cells,
                                                     "T cells" = T_cells,
                                                     "T-rm" = "T-rm",
                                                     "T-NK" = "T-NK"),
                                     features = "IFNG", assays = "SCT"){
        groups = unique(unlist(group.by.list))
        cells_id = (object@meta.data[,group.by] %in% groups) & 
                (object$orig.ident %in% orig.ident)
        cells = colnames(object)[cells_id]
        UMI_df <- FetchData(object = object, vars = c(features,group.by),cells = cells)
        colnames(UMI_df)[1] = "features"
        UMI_list <- list()
        for(i in seq_along(group.by.list)) {
                g = group.by.list[[i]]
                UMI_list[[i]] = UMI_df %>% filter(cell_types %in% g) %>%
                        summarise(all_cells = length(features),
                                  p_positive_cells = sum(features >0)/length(features),
                                  n_positive_cells = sum(features >0),
                                  avg_UMI =log1p(mean(expm1(features))))
                
        }
        UMI =  bind_cols(UMI_list)
        return((UMI))
}
res_list2 <- list()
res_list2[["IFNG immune"]] <- pbapply::pblapply(orig_idents, 
                               FUN = AverageExpression_pct_v2, 
                               features = "IFNG",
                               object=object,
                               group.by.list = list("all immune cells" = Immune_cells,
                                                    "T cells" = T_cells,
                                                    "T-rm" = "T-rm",
                                                    "T-NK" = "T-NK"))
res_list2[["IFNG immune"]] %<>% bind_rows
rownames(res_list2[["IFNG immune"]])  = orig_idents
################## surface airway epithelial ################## 
Idents(SAE) = "cell_types"
jpeg(paste0(path,"IFNGR12_SAE_v1.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(SAE, features = c("IFNGR1","IFNGR2"),group.by = "cell_types", 
        assay = "RNA",pt.size = 0)
dev.off()

# subset
Idents(SAE) = "conditions"
P_D <- subset(SAE,idents = c("proximal","distal"))
P_D$conditions %<>% factor(levels = c("proximal","distal"))
D_COPD <- subset(SAE,idents = c("distal","COPD"))
D_COPD$conditions %<>% factor(levels = c("distal","COPD"))
D <- subset(SAE,idents = "distal")
D$age <- plyr::mapvalues(D$orig.ident, from = age_df$orig.ident, to =age_df$age)
Idents(D) = "age"
D %<>% subset(idents = c("younger","older"))
D$age %<>% factor(levels = c("younger","older"))


jpeg(paste0(path,"IFNGR12_SAE_v2.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(P_D, group.by = "cell_types",features = c("IFNGR1","IFNGR2"),
        assay = "RNA",split.by = "conditions",pt.size = 0)
dev.off()

jpeg(paste0(path,"IFNGR12_SAE_v3.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(D_COPD, group.by = "cell_types",features = c("IFNGR1","IFNGR2"),
        assay = "RNA",split.by = "conditions",pt.size = 0)
dev.off()

jpeg(paste0(path,"IFNGR12_SAE_v4.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(D, group.by = "cell_types",features = c("IFNGR1","IFNGR2"),
        assay = "RNA",split.by = "age",pt.size = 0)
dev.off()

Idents(SAE) = "conditions"
P <- subset(SAE,idents = c("proximal"))
D <- subset(SAE,idents = c("distal"))
COPD <- subset(SAE,idents = c("COPD"))
D$age <- plyr::mapvalues(D$orig.ident, from = age_df$orig.ident, to =age_df$age)
Idents(D) = "age"
younger<- subset(D,idents = "younger")
older<- subset(D,idents = "older")

UMI_list2 <- pbapply::pblapply(list(SAE,P,D,younger,older,COPD), 
                              FUN = AverageExpression_pct,features = "IFNGR1")
UMI_list3 <- pbapply::pblapply(list(SAE,P,D,younger,older,COPD), 
                              FUN = AverageExpression_pct,features = "IFNGR2")
res_list[["IFNGR1 surface airway epi"]] = list2df(UMI_list2)
res_list[["IFNGR2 surface airway epi"]] = list2df(UMI_list3)

openxlsx::write.xlsx(res_list, file =  "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/IFNG-IFNGR data raw.xlsx",
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

########### SAE average expression for per sample ##############
BC_IC = c("BC1", "BC2", "BC-p", "IC1", "IC2", "IC3")
res_list2[["IFNGR1 surface airway epi"]] <- pbapply::pblapply(orig_idents, 
                                                FUN = AverageExpression_pct_v2, 
                                                features = "IFNGR1",
                                                object=object,
                                                group.by.list = list("all SAE_cells" = SAE_cells,
                                                                     "BC_IC" = BC_IC,
                                                                     "BC" = c("BC1", "BC2", "BC-p"),
                                                                     "BC1" = "BC1"))
res_list2[["IFNGR1 surface airway epi"]] %<>% bind_rows
rownames(res_list2[["IFNGR1 surface airway epi"]])  = orig_idents

res_list2[["IFNGR2 surface airway epi"]] <- pbapply::pblapply(orig_idents, 
                                                              FUN = AverageExpression_pct_v2, 
                                                              features = "IFNGR2",
                                                              object=object,
                                                              group.by.list = list("all SAE_cells" = SAE_cells,
                                                                                   "BC_IC" = BC_IC,
                                                                                   "BC" = c("BC1", "BC2", "BC-p"),
                                                                                   "BC1" = "BC1"))
res_list2[["IFNGR2 surface airway epi"]] %<>% bind_rows
rownames(res_list2[["IFNGR2 surface airway epi"]])  = orig_idents

openxlsx::write.xlsx(res_list2, file =  "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/IFNG-IFNGR data per sample raw.xlsx",
                     colNames = TRUE,row.names = T,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

######################### p-value ###########################
deg_list <- pbapply::pblapply(Immune_cells, function(x){
        deg <- readxl::read_excel("Yang/Lung_30/DE_analysis/groups/DE_results_Immune.xlsx",
                                  sheet = x)
        deg[deg$gene %in% "IFNG",]
})
IFNG <- bind_rows(deg_list)

table(Immune$cell_types)%>% prop.table *100 %>% round(digits = 3)

deg_list <- pbapply::pblapply(SAE_cells, function(x){
        deg <- readxl::read_excel("Yang/Lung_30/DE_analysis/groups/DE_results_Surface Airway Epithelial.xlsx",
                                  sheet = x)
        deg[deg$gene %in% c("IFNGR1","IFNGR2"),]
})
IFNGR <- bind_rows(deg_list)
IFNGR$cell_type = gsub(" .*","",IFNGR$cluster)
table(SAE$cell_types) %>% prop.table *100 %>% round(digits = 3)

