library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(magrittr)
# ======== prepare data for Chord diagram, Sankey diagram =============
save.path = "Yang/Lung_30/Cell_Phone_DB/"
# load R-L interaciton
df_res <- readxl::read_excel(paste0(save.path,"Statistics for R-L cell-cell interactions per groups.xlsx"))
colnames.append <- c(rep("",5),rep("_exp",7),rep("_FC",11),rep("_pvalue",11),rep("_p.adj",11))
colnames(df_res) = paste0(df_res[1,], colnames.append)
df_res = df_res[-1,]
df_res1 <- df_res
# rename Cell-cell pair
anno <- readxl::read_excel("doc/Annotations/Cell type abbreviation.xlsx")
cell.types = anno$Abbreviation
cell_types = gsub('-','_',cell.types)
for(i in grep("-",cell.types)){
    df_res$`Cell-cell pair` %<>% gsub(paste0("^",cell.types[i],"-"),paste0(cell_types[i],"-"),.)
    df_res$`Cell-cell pair` %<>% gsub(paste0("-",cell.types[i],"$"),paste0("-",cell_types[i]),.)
    svMisc::progress(i,length(cell.types))
}
df_res$`Cell-cell pair` %<>% gsub(paste0("_p-C$"),"-p_C",.)
df_res$`Cell-cell pair` %<>% gsub(paste0("_S-d$"),"-S_d",.)

df_res %<>% cbind(stringr::str_split(string = df_res$`Cell-cell pair`,pattern = "-",simplify = T))
df_res$`1` %<>% gsub("_","-",.)
df_res$`2` %<>% gsub("_","-",.)
df_res$`1` %<>% plyr::mapvalues(from = anno$Abbreviation,
                                to = anno$`Revised abbreviations`)
df_res$`2` %<>% plyr::mapvalues(from = anno$Abbreviation,
                                to = anno$`Revised abbreviations`)
df_res$`Cell-cell pair` = paste0(df_res$`1`,"_", df_res$`2`)
# subset R-L interaciton data

select_col = c("R","L","R-L pair","Cell-cell pair","P_exp","P vs D+T_FC","P vs D+T_pvalue")
sub_res = df_res[,select_col]
sub_res[,select_col[5:7]] %<>% lapply(function(x) as.numeric(x))
# filtering
select_row = sub_res[,"P vs D+T_FC"] > 1 & sub_res[,"P vs D+T_pvalue"] < 0.05
table(select_row)
sub = sub_res[select_row,]
sub = data.table::data.table(sub)
sub = sub[,c("R-L pair","Cell-cell pair")]
sub = sub[,list(count=.N),by=`Cell-cell pair`]
dim(sub)
# cell order  
df <- readxl::read_excel("doc/Chord diagram cell order.xlsx",col_names = T)
select_celltypes = df$cell.types
anno$`Revised abbreviations`[!(anno$`Revised abbreviations` %in% select_celltypes)]

# prepare cell-cell pair with 0 count
All_cell_pair = paste0(rep(df$cell.types, each = 59), "_",rep(df$cell.types, time = 59))
count_0 = F
if(count_0){
    table(All_cell_pair %in% sub$`Cell-cell pair`)
    sub1 = data.frame(`Cell-cell pair` = All_cell_pair[!(All_cell_pair %in% sub$`Cell-cell pair`)],
                           `count` = 0)
    colnames(sub1)[1] = "Cell-cell pair"
    sub %<>% as.data.frame %>% rbind(sub1)
}


sub %<>% as.data.frame %>% cbind(str_split(string = sub$`Cell-cell pair`,pattern = "_",simplify = T))
table(keep <- (sub$`1` %in% select_celltypes) & (sub$`2` %in% select_celltypes))
sub = sub[keep,]
sub$`1` %<>% gdata::reorder.factor(new.order=select_celltypes)
sub %<>%   arrange(`1`)
table(sub$`Cell-cell pair` %in% All_cell_pair)

#  ============== Alluvial Plots =====================
library(ggalluvial)

# I need a long format
sub_long = gather(sub, key = "signal_receiver", value = "Cell.types",
                  -c("Cell-cell pair","count"))
grid.col = df$hexcolor[select_celltypes %in% unique(pull(sub[,c("1","2")]))]

sub_long$signal_receiver %<>% plyr::mapvalues(from = c("1", "2"),
                                              to = c("Ligands","Receptors"))
sub_long$signal_receiver %<>% as.factor()
sub_long$signal_receiver %<>% factor(levels = c("Ligands","Receptors"))

sub_long$Cell.group =  plyr::mapvalues(x = sub_long$Cell.types,
                                       from = df$cell.types,
                                       to = df$Cell.group)
sub_long$Cell.types %<>% gdata::reorder.factor(new.order=select_celltypes)
sub_long %<>%   arrange(Cell.types)
sub_long$Cell.group %<>% as.factor()
sub_long$Cell.group %<>% factor(levels = c("Epithelial","Structural","Immune"))

g <- ggplot(sub_long,
            aes(x = signal_receiver, stratum = Cell.types, alluvium = `Cell-cell pair`,
                y = count,
                fill = Cell.group, label = Cell.types)) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    #geom_text(stat = "stratum", size = 3,nudge_x = 0.25) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.6)) +
    scale_fill_manual(values=unique(df$hexcolor))+
    ggtitle(paste("Ligands-Receptor Signal in","P vs D+T"))

pbuild <- ggplot_build(plot = g)
x.min <- min(pbuild$layout$panel_params[[1]]$x.range) + 0.1
x.max <- max(pbuild$layout$panel_params[[1]]$x.range) - 0.1
g <- g +     coord_cartesian(xlim = c(x.min, x.max), clip = "off") 

# Make the Network
jpeg(paste0(save.path,"alluvium_","P vs D+T_nolabel",".jpeg"), family = "Arial",
     units="in", width=7, height=25,res=900)
print(g)
dev.off()

##==========  Chord diagram ==================
# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggraph)
library(colormap)
library(tibble)
library(igraph)
library(spatstat)

# Transform the adjacency matrix in a long format
sub1 <-  sub
sub <- sub1
colnames(sub) = c("pair","value","from","to")
sub = sub[,c("from","to","value","pair")]
# Number of connection per person
    as_tibble(sub) %>%
    group_by(to) %>%
    summarize(sum=sum(value)) -> coauth
colnames(coauth) <- c("name", "sum")
table(coauth$sum)
hist(coauth$sum)
# Create a graph object with igraph
mygraph <- graph_from_data_frame(d = sub, vertices = coauth, directed = TRUE )

# Find community
#com <- walktrap.community(mygraph)
com <- df[match(coauth$name, df$cell.types),"Cell.group"]
com = na.omit(com$Cell.group)
#Reorder dataset and make the graph

coauth <- coauth %>% 
    mutate( grp = com) %>%
    arrange(grp) %>%
    mutate(name=factor(name, name))

# keep only 10 first communities
#coauth <- coauth %>% 
#    filter(sum > 75)

# keep only this people in edges
sub <- sub %>%
    filter(from %in% coauth$name) %>%
    filter(to %in% coauth$name)

table(c(as.character(sub$from), as.character(sub$to)) %in% coauth$name)
table(coauth$name %in% c(as.character(sub$from), as.character(sub$to)))

# Add label angle
number_of_bar=nrow(coauth)
coauth$id = seq(1, nrow(coauth))
angle= 360 * (coauth$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
coauth$hjust <- ifelse(angle > 90 & angle<270, 1, 0)
coauth$angle <- ifelse(angle > 90 & angle<270, angle+180, angle)

# Create a graph object with igraph
mygraph <- graph_from_data_frame( sub, vertices = coauth, directed = TRUE )

# prepare a vector of n color in the viridis scale
mycolor <- unique(df$hexcolor)[match(unique(df$Cell.group),unique(coauth$grp))]
# Make the graph

g <- ggraph(mygraph, layout="circle") + 
    geom_edge_link(edge_colour="black", edge_alpha=0.2, edge_width=0.3) +
    geom_node_point(aes(size=sum, color=as.factor(grp), fill=grp), alpha=0.9) +
    scale_size_continuous(range=c(0.5,8)) +
    scale_color_manual(values=mycolor) +
    geom_node_text(aes(label=paste("    ",name,"    "), angle=angle, hjust=hjust), size=2.3, color="black") +
    theme_void() +
    theme(
        legend.position="none",
        plot.margin=unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) 


jpeg(paste0(save.path,"new_chordDiagram_","P vs D+T",".jpeg"), family = "Arial",
     units="in", width=7, height=7,res=900)
print(g)
dev.off()
#========== circlize =============

# https://www.data-to-viz.com/story/AdjacencyMatrix.html
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag)  #devtools::install_github("mattflor/chorddiag")
library(magrittr)
library(data.table)
# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
mycolor <- viridis(10, alpha = 1, begin = 0, end = 1, option = "D")
mycolor <- mycolor[sample(1:10)]

# Base plot
jpeg(paste0(save.path,"chordDiagram_","P vs D+T",".jpeg"), family = "Arial",
     units="in", width=7, height=7,res=900)
chordDiagram(
    x = as.data.frame(sub[,c("V1","V2","count")]), 
    grid.col = df$hexcolor[select_celltypes %in% unique(pull(sub[,c("V1","V2")]))],
    transparency = 0.25,
    directional = 1,
    direction.type = c("arrows", "diffHeight"), 
    diffHeight  = -0.04,
    #annotationTrack = "grid", 
    annotationTrackHeight = c(0.05, 0.1),
    link.arr.type = "big.arrow", 
    link.sort = F, 
    link.largest.ontop = F)
dev.off()

# Add text and axis
circos.trackPlotRegion(
    track.index = 1, 
    bg.border = NA, 
    panel.fun = function(x, y) {
        
        xlim = get.cell.meta.data("xlim")
        sector.index = get.cell.meta.data("sector.index")
        
        # Add names to the sector. 
        circos.text(
            x = mean(xlim), 
            y = 3.2, 
            labels = sector.index, 
            facing = "bending", 
            cex = 0.8
        )
        
        # Add graduation on axis
        circos.axis(
            h = "top", 
            major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
            minor.ticks = 1, 
            major.tick.length = 0.5,
            labels.niceFacing = FALSE)
    }
)

#  ============== Sankey diagram =====================
library(ggalluvial)
library(networkD3)
# I need a long format
sub_long = as.data.frame(sub[,c("V1","V2","count")])
colnames(sub_long) <- c("source", "target", "value")
sub_long$target <- paste(sub_long$target, " ", sep="")

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(sub_long$source), as.character(sub_long$target)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
sub_long$IDsource=match(sub_long$source, nodes$name)-1 
sub_long$IDtarget=match(sub_long$target, nodes$name)-1

# prepare colour scale
grid.col = df$hexcolor[select_celltypes %in% unique(pull(sub[,c("V1","V2")]))]
grid.col = paste0("d3.scaleOrdinal() .range(['",paste(grid.col,collapse = "','"),"'])")

ColourScal ="d3.scaleOrdinal() .range(['#FDE725FF','#B4DE2CFF','#6DCD59FF','#35B779FF','#1F9E89FF','#26828EFF','#31688EFF','#3E4A89FF','#482878FF','#440154FF'])"

# Make the Network
jpeg(paste0(save.path,"sankeyNetwork_","P vs D+T",".jpeg"), family = "Arial",
     units="in", width=7, height=21,res=900)
sankeyNetwork(Links = sub_long, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=T, #colourScale=grid.col,
              nodeWidth=40, fontSize=13, nodePadding=20)
dev.off()
