library(ggpubr)
library(tidyr)
library(magrittr)
library(stringr)
library(ggExtra)
library(patchwork)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

data <- readxl::read_excel("Yang/Lung_30/Monocle2/Pseudotime plot-1.xlsx", sheet = "Sheet1")
sample_df <- data[,c("Regions","Samples")]
df <- data[complete.cases(data),-c(1,2)]
df$Pseudotime %<>% paste0("_",df$State)
df <- gather(data = df[,-1],key ="samples",value = "value", -Pseudotime)
df$log_value = log1p(df$value)
Pseudotime_state <- str_split(df$Pseudotime, pattern = "_")
df$Pseudotime = sapply(Pseudotime_state,function(x) as.integer(x[[1]]))
df$State = sapply(Pseudotime_state,function(x) x[[2]]) %>% 
        factor(levels = paste("State",c(1,6,7,15,16,17)))
df$Regions = plyr::mapvalues(df$samples, from = sample_df$Samples, to = sample_df$Regions)
df$Regions %<>% plyr::mapvalues(from = c("P","D","T"),
                                to = c("Proximal","Distal","Terminal"))
df$Regions %<>% factor(levels = c("Proximal","Distal","Terminal","COPD"))

jpeg(paste0(path,"Pseudotime_log_PDT.jpeg"), units="in", width=10, height=7,res=600)
ggscatter(df[df$Regions %in% c("Proximal","Distal","Terminal"),],
          x = "Pseudotime", y = "value",
          color = "Regions",palette = c("#1F78B4","#4ca64c","#E6AB02"),
          shape = "State",                        # Extending the regression line
          add = "loess",conf.int = TRUE)+
        stat_cor(aes(color = Regions), method = "spearman", label.x = 3)           # Add correlation coefficient

dev.off()

jpeg(paste0(path,"Pseudotime_log_CODP.jpeg"), units="in", width=10, height=7,res=600)
ggscatter(df[df$Regions %in% c("Distal","COPD"),],
          x = "Pseudotime", y = "log_value",
          color = "Regions",palette = c("#4ca64c","#C53B19"),
          shape = "State",                        # Extending the regression line
          add = "loess",conf.int = TRUE)+
        stat_cor(aes(color = Regions), method = "spearman", label.x = 3)           # Add correlation coefficient
dev.off()

jpeg(paste0(path,"Pseudotime_PDT.jpeg"), units="in", width=10, height=7,res=600)
ggscatter(df[df$Regions %in% c("Proximal","Distal","Terminal"),],
          x = "Pseudotime", y = "value",
          color = "Regions",palette = c("#1F78B4","#4ca64c","#E6AB02"),
          shape = "State",                        # Extending the regression line
          add = "loess",conf.int = TRUE)+
        stat_cor(aes(color = Regions), method = "spearman", label.x = 3)           # Add correlation coefficient
dev.off()

# ==================== ggridges ============
# try to construct a fabricated data to produce Density ridgeline plots
df$value %<>% ceiling()
#mtx$barcode = stringi::stri_rand_strings(n = nrow(mtx), length = 6, '[A-Z]')

mtx <- df[rep(1, each = as.integer(df[1,"value"])), ]
for(i in 2:nrow(df)) {
        mtx %<>% rbind(df[rep(i, each = as.integer(df[i,"value"])), ])
        svMisc::progress(i/nrow(df)*100)
}

g_list <- pbapply::pblapply(c("Proximal","Distal","Terminal"), function(x) {
        g <- ggplot(mtx[mtx$Regions %in% x,], 
                    aes(x = Pseudotime,y = State, fill = stat(x)))+ 
                geom_density_ridges_gradient(rel_min_height = 0.01, scale = 3, alpha = 0.7)+
                scale_x_continuous(expand = c(0, -1)) +
                scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
                scale_fill_viridis_c(name = "Pseudotime", option = "D")+
                theme_ridges()+
                labs(title = x)+ 
                TitleCenter(size = 20)+
                NoLegend()
        keep.text <- x == "Proximal"
        g + NoAxes(keep.text = keep.text)
})

jpeg(paste0(path,"Pseudotime_PDT.jpeg"), units="in", width=10, height=7,res=600)
wrap_plots(g_list, nrow = 1)
dev.off()
#Blue (P): 
#Filling: 137, 207, 240 #89CFF0
#Line: 0, 125, 220 #007DDC

#Green (D): 
#Filling: 162, 241, 155 #A2F19B
#Line: 31, 198, 80 #1FC650

#Yellow (T):
#Filling: 255, 244, 153 #FFF499
#Line:  255, 206, 0 #FFCE00

#Red (COPD):
#Filling: 254, 200, 200 #FEC8C8
#Line: 239, 71, 70 #EF4746
g_list <- list()
for(i in seq_along(c("Proximal","Distal","Terminal"))) {
        x = c("Proximal","Distal","Terminal")[i]
        g <- ggplot(mtx[mtx$Regions %in% x,], 
                    aes(x = Pseudotime,y = State, fill = x, color = x))+ 
                geom_density_ridges_gradient(rel_min_height = 0.01, scale = 3, 
                                             alpha = 0.7)+
                scale_fill_cyclical(values = c("#89CFF0","#A2F19B","#FFF499")[i], guide = "none") +
                scale_color_cyclical(values = c("#007DDC","#1FC650","#FFCE00")[i],guide = "none") +
                scale_x_continuous(expand = c(0, -1)) +
                scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
                theme_ridges()+
                labs(title = x)+ 
                TitleCenter(size = 20)+
                NoLegend()
        keep.text <- x == "Proximal"
        g_list[[i]] = g + NoAxes(keep.text = keep.text)
}
jpeg(paste0(path,"Pseudotime_PDT~.jpeg"), units="in", width=10, height=7,res=600)
wrap_plots(g_list, nrow = 1)
dev.off()


g_list <- pbapply::pblapply(c("Distal","COPD"), function(x) {
        g <- ggplot(mtx[mtx$Regions %in% x,], 
                    aes(x = Pseudotime,y = State, fill = x, color = x))+ 
                geom_density_ridges_gradient(rel_min_height = 0.01, scale = 3, alpha = 0.7)+
                scale_x_continuous(expand = c(0, -1)) +
                scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
                scale_fill_viridis_c(name = "Pseudotime", option = "D")+
                theme_ridges()+
                labs(title = x)+ 
                TitleCenter(size = 20)+
                NoLegend()
        keep.text <- x == "Distal"
        g + NoAxes(keep.text = keep.text)
})

jpeg(paste0(path,"Pseudotime_CODP.jpeg"), units="in", width=10, height=7,res=600)
wrap_plots(g_list, nrow = 1)
dev.off()
#Green (D): 
#Filling: 162, 241, 155 #A2F19B
#Line: 31, 198, 80 #1FC650

#Blue (P): 
#Filling: 137, 207, 240 #89CFF0
#Line: 0, 125, 220 #007DDC

#Yellow (T):
#Filling: 255, 244, 153 #FFF499
#Line:  255, 206, 0 #FFCE00

#Red (COPD):
#Filling: 254, 200, 200 #FEC8C8
#Line: 239, 71, 70 #EF4746

g_list <- list()
for(i in seq_along(c("Distal","COPD"))) {
        x = c("Distal","COPD")[i]
        g <- ggplot(mtx[mtx$Regions %in% x,], 
                    aes(x = Pseudotime,y = State, fill = x, color = x))+ 
                geom_density_ridges_gradient(rel_min_height = 0.01, scale = 3, alpha = 0.7)+
                scale_fill_cyclical(values = c("#A2F19B","#FEC8C8")[i], guide = "none") +
                scale_color_cyclical(values = c("#1FC650","#EF4746")[i],guide = "none") +
                scale_x_continuous(expand = c(0, -1)) +
                scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
                theme_ridges()+
                labs(title = x)+ 
                TitleCenter(size = 20)+
                NoLegend()
        keep.text <- x == "Proximal"
        g_list[[i]] = g + NoAxes(keep.text = keep.text)
}
jpeg(paste0(path,"Pseudotime_CODP~.jpeg"), units="in", width=10, height=7,res=600)
wrap_plots(g_list, nrow = 1)
dev.off()


g <- ggplot(mtx[mtx$Regions %in% c("Proximal","Distal","Terminal"),],
            aes(x = Pseudotime, y = State, color = Regions, fill = NA)) +
        geom_density_ridges(scale = 2, rel_min_height = .01,alpha = 0.7,size=1)+
        scale_x_continuous(expand = c(0, -1)) +
        scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
        scale_fill_manual(values = c("#007DDC","#1FC650", "#FFCE00"),
                          labels = c("Proximal","Distal","Terminal")) +
        scale_color_manual(values = c("#007DDC","#1FC650",  "#FFCE00"), guide = "none") +
        ggtitle("Compare pseudotime in Proximal,Distal, and Terminal") +
        theme_ridges(center = TRUE)+ TitleCenter()

jpeg(paste0(path,"Pseudotime_PDT_combine.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()


g <- ggplot(mtx[mtx$Regions %in% c("Distal","COPD"),],
            aes(x = Pseudotime, y = State, color = Regions, fill = NA)) +
        geom_density_ridges(scale = 2, rel_min_height = .01,alpha = 0.7,size=1)+
        scale_x_continuous(expand = c(0, -1)) +
        scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
        scale_fill_manual(values = c("#1FC650", "#EF4746"),
                          labels = c("Distal","COPD")) +
        scale_color_manual(values = c("#1FC650", "#EF4746"), guide = "none") +
        ggtitle("Compare pseudotime in Proximal,Distal, and Terminal") +
        theme_ridges(center = TRUE)+ TitleCenter()

jpeg(paste0(path,"Pseudotime_CODP_combine.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()
