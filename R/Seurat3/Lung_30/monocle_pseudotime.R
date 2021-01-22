library(ggpubr)
library(tidyr)
library(magrittr)
library(stringr)
library(ggExtra)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

data <- readxl::read_excel("Yang/Lung_30/Monocle2/Pseudotime plot-1.xlsx", sheet = "Sheet1")
sample_df <- data[,c("Regions","Samples")]
df <- data[complete.cases(data),-c(1,2)]
df$Pseudotime %<>% paste0("_",df$State)
df <- gather(data = df[,-1],key ="samples",,  -Pseudotime)
df$log_value = log1p(df$value)
Pseudotime_state <- str_split(df$Pseudotime, pattern = "_")
df$Pseudotime = sapply(Pseudotime_state,function(x) as.integer(x[[1]]))
df$State = sapply(Pseudotime_state,function(x) x[[2]])
df$Regions = plyr::mapvalues(df$samples, from = sample_df$Samples, to = sample_df$Regions)
df$Regions %<>% factor(levels = c("P","D","T","COPD"))
df1 <- df[df$Regions %in% c("P","D","T"),]
df2  <- df[df$Regions %in% c("D","COPD"),]

jpeg(paste0(path,"Pseudotime_log_PDT.jpeg"), units="in", width=10, height=7,res=600)
ggscatter(df1, x = "Pseudotime", y = "log_value",
          color = "Regions",palette = c("#1F78B4","#4ca64c","#E6AB02"),
          shape = "State",                        # Extending the regression line
          add = "loess",conf.int = TRUE)+
        stat_cor(aes(color = Regions), label.x = 3)           # Add correlation coefficient
dev.off()

jpeg(paste0(path,"Pseudotime_log_CODP.jpeg"), units="in", width=10, height=7,res=600)
ggscatter(df2, x = "Pseudotime", y = "log_value",
          color = "Regions",palette = c("#4ca64c","#C53B19"),
          shape = "State",                        # Extending the regression line
          add = "loess",conf.int = TRUE)+
        stat_cor(aes(color = Regions), method = "spearman", label.x = 3)           # Add correlation coefficient
dev.off()

jpeg(paste0(path,"Pseudotime_PDT.jpeg"), units="in", width=10, height=7,res=600)
ggscatter(df1, x = "Pseudotime", y = "value",
          color = "Regions",palette = c("#1F78B4","#4ca64c","#E6AB02"),
          shape = "State",                        # Extending the regression line
          add = "loess",conf.int = TRUE)+
        stat_cor(aes(color = Regions), method = "spearman", label.x = 3)           # Add correlation coefficient
dev.off()
