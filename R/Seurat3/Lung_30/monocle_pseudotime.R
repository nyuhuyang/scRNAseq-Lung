library(ggpubr)
library(tidyr)
library(magrittr)
library(stringr)
library(ggExtra)
library(patchwork)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

data <- readxl::read_excel("Yang/Lung_30/Monocle2/Pseudotime plot-1.xlsx", sheet = "Sheet1")
sample_df <- data[,c("Regions","Samples")]
df <- data[complete.cases(data),-c(1,2)]

#df$Pseudotime_state = paste0(df$Pseudotime,"_",df$State)
df <- gather(data = df,key ="samples",value = "value", -c("Pseudotime","State","Pseudotime1"))
df$log_value = log1p(df$value)
#Pseudotime_state <- str_split(df$Pseudotime_state, pattern = "_")
#df$Pseudotime1 = df$Pseudotime
#df$Pseudotime = sapply(Pseudotime_state,function(x) as.integer(x[[1]]))
#= sapply(Pseudotime_state,function(x) x[[2]]) 
df$State %<>% factor(levels = paste("State",c(15,7,6,16,17,1)))
df$Regions = plyr::mapvalues(df$samples, from = sample_df$Samples, to = sample_df$Regions)
df$Regions %<>% plyr::mapvalues(from = c("P","D","T"),
                                to = c("Proximal","Distal","Terminal"))
df$Regions %<>% factor(levels = c("Proximal","Distal","Terminal","COPD"))
df$Pseudotime1 %<>% as.integer()

# ========  ggbarplot   ===========
df1 <- df %>% filter(Regions  %in% c("Proximal","Terminal","Distal"))
g <- ggboxplot(data = df1,
          x = "Pseudotime1", y ="value",
          #add = c("mean_sd"),
          bxp.errorbar = T,
          bxp.errorbar.width = 0.4,
          error.plot = "errorbar",
          width=0.7,
          color = "Regions", 
          fill = "Regions",
          outlier.shape = 21,
          #position = position_dodge(0.75),
          title = "Compare pseudotime in Proximal, Distal, and Terminal",
          xlab = "Pseudotime_code",
          ylab ="cell percentage (%)",
          palette = c("#007DDC","#1FC650",  "#FFCE00"))+
        scale_fill_manual(values = c("#89CFF0","#A2F19B", "#FFF499"),
                          labels = c("Proximal","Distal","Terminal"))+
        geom_smooth(aes(group=Regions, color = Regions),
                    method = 'loess', formula = y ~ x, se=FALSE,
                    size = 0.5)+ TitleCenter()

        
jpeg(paste0(path,"Pseudotime_PDT_dots.jpeg"), units="in", width=12, height=7,res=600)
print(g)
dev.off()


df2 <- df %>% filter(Regions  %in% c("Distal","COPD"))
g <- ggboxplot(data = df2,
               x = "Pseudotime1", y ="value",
               #add = c("mean_sd", "point"),
               bxp.errorbar = T,
               bxp.errorbar.width = 0.4,
               error.plot = "errorbar",
               width=0.7,
               color = "Regions", 
               fill = "Regions",
               outlier.shape = 21,
               #position = position_dodge(0.75),
               title = "Compare pseudotime in Distal and COPD",
               xlab = "Pseudotime_code",
               ylab ="cell percentage (%)",
               palette = c("#1FC650", "#EF4746"))+
        scale_fill_manual(values = c("#A2F19B", "#FEC8C8"),
                          labels = c("Distal","COPD"))+
        geom_smooth(aes(group=Regions, color = Regions),
                    method = 'loess', formula = y ~ x, se=FALSE,
                    size = 0.5)+ TitleCenter()


jpeg(paste0(path,"Pseudotime_COPD_dots.jpeg"), units="in", width=12, height=7,res=600)
print(g)
dev.off()
#====

g <- ggplot(data = df %>% filter(Regions  %in% c("Proximal","Terminal","Distal")) %>%
                    arrange(desc(Pseudotime1)),
        aes(x = Pseudotime1, y =value))
g + geom_boxplot(aes(fill = Regions), width = 0.5, size = 0.4,
                     position = position_dodge(0.7),
                 fill=Regions, color = Regions,
               x = "Pseudotime1", y ="value",
               add = c("mean_sd", "point"),
               error.plot = "errorbar",
               add.params = list(width=1, size = 0.5),

               position = position_dodge(0.7),
               palette = c("#007DDC","#1FC650",  "#FFCE00"))+
        scale_fill_manual(values = c("#89CFF0","#A2F19B", "#FFF499"),
                          labels = c("Proximal","Distal","Terminal"))+
        geom_smooth(data = df %>% filter(Regions  %in% c("Proximal","Terminal","Distal")),
                    aes(group=Regions, color = Regions),
                    method = 'loess', formula = y ~ x,se=FALSE,
                    size = 0.5)+
        ggtitle("Compare pseudotime in Proximal, Distal, and Terminal") +
        TitleCenter()+
        xlab("Pseudotime_code")+
        ylab("cell percentage (%)")


jpeg(paste0(path,"Pseudotime_PDT_dots.jpeg"), units="in", width=15, height=7,res=600)
print(g)
dev.off()



# ========  ggscatter   ===========
jpeg(paste0(path,"Pseudotime_PDT.jpeg"), units="in", width=10, height=7,res=600)
ggscatter(df[df$Regions %in% c("Proximal","Distal","Terminal"),],
          x = "Pseudotime1", y = "log_value",
          color = "Regions",
          shape = NA,                        # Extending the regression line
          palette = c("#007DDC","#1FC650",  "#FFCE00"),   
          add = "loess",add.params =  list(size=0.5),
          conf.int = TRUE)+
        stat_cor(aes(color = Regions), method = "spearman", label.x = 3)+           # Add correlation coefficient
        scale_x_continuous(name ="Pseudotime",
                           breaks = df$Pseudotime1,
                           labels = df$Pseudotime)
dev.off()

jpeg(paste0(path,"Pseudotime_CODP.jpeg"), units="in", width=10, height=7,res=600)
ggscatter(df[df$Regions %in% c("Distal","COPD"),],
          x = "Pseudotime1", y = "log_value",
          color = "Regions",palette = c("#1FC650", "#EF4746"),
          shape = NA,                        # Extending the regression line
          add = "loess",add.params = list(size=0.5),
          conf.int = TRUE)+
        stat_cor(aes(color = Regions), method = "spearman", label.x = 3)+           # Add correlation coefficient
        scale_x_continuous(name ="Pseudotime",
                           breaks = df$Pseudotime1,
                           labels = df$Pseudotime)
dev.off()

# ==================== ggridges ============
# try to construct a fabricated data to produce Density ridgeline plots
df1 <- df %>% group_by(Regions, State, Pseudotime) %>% 
        summarise(mean_value =mean(value)) %>%
                ungroup()
#df1$value = log1p(df1$mean_value *1)
#df1[df1$value == 0,"value"]  = df1[df1$value == 0,"value"] +1
#df1 <- df

mtx_list <- list()
for(i in 1:nrow(df1)) {
        mtx_list[[i]] = df1[rep(i, each = as.integer(df1[i,"mean_value"])), ]
        svMisc::progress(i/nrow(df1)*100)
}
mtx = bind_rows(mtx_list)

g <- mtx %>% filter(Regions  %in% c("Proximal","Terminal","Distal")) %>%
        gghistogram(x = "Pseudotime",
                    #add = "mean", #rug = TRUE,
                    position = "dodge",
                    fill = "Regions", palette = c("#007DDC","#1FC650",  "#FFCE00"),
                    add_density = T)+
        scale_fill_manual(values = c("#89CFF0","#A2F19B", "#FFF499"),
                          labels = c("Proximal","Distal","Terminal"))+
        ggtitle("Compare pseudotime in Proximal,Distal, and Terminal") +
        TitleCenter()+
        ylab("Value")


jpeg(paste0(path,"Pseudotime_PDT.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()

g <- mtx %>% filter(Regions  %in% c("Distal","COPD")) %>%
        gghistogram(x = "Pseudotime",
                    #add = "mean", #rug = TRUE,
                    position = "dodge",
                    fill = "Regions", palette = c("#1FC650", "#EF4746"),
                    add_density = TRUE)+
        scale_fill_manual(values = c("#A2F19B", "#FEC8C8"),
                          labels = c("Distal","COPD"))+
        ylab("Value")+
        ggtitle("Compare pseudotime in Distal and COPD") +
        TitleCenter()
                    
jpeg(paste0(path,"Pseudotime_CODP.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()
#Green (D): 
#Filling: 162, 241, 155 #A2F19B
#Line: 31, 198, 80 #1FC650

#Red (COPD):
#Filling: 254, 200, 200 #FEC8C8
#Line: 239, 71, 70 #EF4746





library(ggridges)
g_list <- pbapply::pblapply(c("Proximal","Distal","Terminal"), function(x) {
        g <- ggplot(mtx[mtx$Regions %in% x,], 
                    aes(x = Pseudotime,y = State, fill = stat(x)))+ 
                geom_density_ridges_gradient(rel_min_height = 0.01, scale = 1, alpha = 0.7)+
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
                geom_density_ridges_gradient(rel_min_height = 0.01, scale = 1, 
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
                theme(text = element_text(size=20),
                      plot.title = element_text(hjust = 0.5))+
                NoLegend()
        keep.text <- x == "Proximal"
        g_list[[i]] = g + NoAxes(keep.text = keep.text)
}
jpeg(paste0(path,"Pseudotime_CODP~.jpeg"), units="in", width=10, height=7,res=600)
wrap_plots(g_list, nrow = 1)
dev.off()


mtx %>% filter(Regions  %in% c("Proximal","Terminal","Distal")) %>%
        filter(State %in% c("State 17","State 1")) %>%
        ggplot()+
        aes(x = Pseudotime, y = Regions, height = ..density.., color = Regions, fill = Regions) +
        geom_density_ridges(stat="density_ridges", scale = 1, rel_min_height = .01,alpha = 0.3,size=1)+
        scale_x_continuous(expand = c(0, -1)) +
        scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
        #scale_fill_manual(values = c("#89CFF0","#A2F19B", "#FFF499"),
        #                  labels = c("Proximal","Distal","Terminal")) +
        #scale_color_manual(values = c("#007DDC","#1FC650",  "#FFCE00"), guide = "legend") +
        ggtitle("Compare pseudotime in Proximal,and Terminal") +
        theme_ridges(center = TRUE)+ TitleCenter()
print(g)

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

#================= test ====================
mtx %>% filter(Regions  %in% c("Proximal","Terminal","Distal")) %>%
        filter(State %in% c("State 17","State 1")) %>%
        #filter(samples %in% c("UNC-48-P","UNC-48","CU-12")) %>%
        ggscatter(x = "Pseudotime", y = "Regions",height = "Regions",
                  color = "Regions",
                  shape = NA,                        # Extending the regression line
                  palette = c("#007DDC","#1FC650",  "#FFCE00"))+
        #aes(x = Pseudotime, y = value, height = value, color = Regions, fill = Regions) +
        geom_density_ridges(stat="density_ridges", scale = 1, 
                            rel_min_height = 0,alpha = 0.3,size=1)+
        scale_x_continuous(expand = c(0, -1)) +
        #scale_y_discrete(expand = expansion(mult = c(0, .7))) +

        #scale_color_manual(values = c("#007DDC","#1FC650",  "#FFCE00"), guide = "legend") +
        ggtitle("Compare pseudotime in Proximal, Distal, and Terminal") +
        theme_ridges(center = TRUE)+ TitleCenter()


mtx %>% filter(Regions  %in% c("Proximal","Terminal","Distal")) %>%
        #filter(State %in% c("State 17","State 1")) %>%
        gghistogram(x = "Pseudotime",
                    add = "mean", #rug = TRUE,
                    position = "dodge",
                    fill = "Regions", palette = c("#007DDC","#1FC650",  "#FFCE00"),
                    add_density = TRUE)+
        scale_fill_manual(values = c("#89CFF0","#A2F19B", "#FFF499"),
                          labels = c("Proximal","Distal","Terminal"))+
        ylab("Value")




