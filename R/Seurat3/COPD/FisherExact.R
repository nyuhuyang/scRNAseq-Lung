library(readxl)
library(dplyr)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

(load(file = "data/Lung_16_distal_20191022.Rda"))

# Fisher's Exact Test for values in Seurat
FisherTestforSeurat <- function(object,var = c("singler1main","orig.ident"),percentage = FALSE,
                             file.name = "cell_numbers",...){
        (df <- table(object@meta.data[,var[1]],object@meta.data[,var[2]]) %>% as.data.frame %>%
                 spread(Var2, Freq))
        colnames(df)[1] = "Cell.type"
        rownames(df) = df$Cell.type
        df = df[,-1]
        
        if(percentage) {
                (df1 <- table(object@meta.data[,var[1]],object@meta.data[,var[2]]) %>% 
                         prop.table(2) %>%
                         as.data.frame %>%
                         spread(Var2, Freq))
                colnames(df1)[1] = "Cell.type"
                rownames(df1) = df1$Cell.type
                df1 = df1[,-1]
        }        
        p_value = c()
        ColSum <- colSums(df)
        
        for(i in 1:nrow(df)){
                (conting <- rbind(df[i,],ColSum-df[i,]))
                FISH <- fisher.test(conting,workspace = 2000000,...)
                (p_value[i] = FISH$p.value)
                #CHI = chisq.test(conting, correct = T)
                #chisq_p_value[i] = CHI$p.value             
        }
        
        df$p_value = p_value
        df$p_val_adj = p.adjust(p = df$p_value, method = "bonferroni", 
                                n = nrow(df))
        if(percentage) df = cbind(df1,df)
        write.csv(df,paste0(path,file.name,".csv"))
        return(df)
}
FisherTestforSeurat(object,var = c("labels","conditions"),
                 percentage= F, file.name = "cell_numbers") %>% kable %>% kable_styling

FisherTestforSeurat(object,var = c("labels","conditions"),
                 percentage= T, file.name = "cell_numbers_p") %>% 
        kable %>% kable_styling
