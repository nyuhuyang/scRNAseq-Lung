library(stringr)
library(dplyr)
library(tibble)
library(magrittr)
# Chord diagram  Libraries
library(hrbrthemes)
library(viridis)
library(ggraph)
library(colormap)
library(tibble)
library(igraph)
library(spatstat)
library(RColorBrewer)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

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

# cell order and hierarchy
df <- readxl::read_excel("doc/Chord diagram cell order - updated abbreviations 12-14-20.xlsx",col_names = T)
select_celltypes = df$cell.types
anno$`Revised abbreviations`[!(anno$`Revised abbreviations` %in% select_celltypes)]

cell_order_list <- df %>% split(df$Cell.group)
cell_order = c(rev(cell_order_list[["Epithelial"]]$cell.types)[8:22],
               rev(cell_order_list[["Immune"]]$cell.types),
               rev(cell_order_list[["Structural"]]$cell.types),
               rev(cell_order_list[["Epithelial"]]$cell.types)[1:7])
# Origin on top, then section, then subgroups
d1 <- data.frame(Cell.group="origin", cell.types=c("Epithelial","Structural","Immune"))
hierarchy <- rbind(d1, df[,c("Cell.group","cell.types")])
# create a vertices data.frame. One line per object of our hierarchy, giving features of nodes.
vertices <- data.frame(name = unique(c(as.character(hierarchy$Cell.group), as.character(hierarchy$cell.types))) ) 
#
d = df_res
select ="P vs D+T"
FC = 4
count_0 = T
Hierarchical_Edge_Bundling <- function(d = df_res, select ="P vs D+T", FC = 2, count_0 = T,save.path= NULL){
    # https://www.r-graph-gallery.com/311-add-labels-to-hierarchical-edge-bundling.html
    # subset R-L interaciton data
    select_col = c("R","L","R-L pair","Cell-cell pair",
                   paste0(sub(" .*","",select),"_exp"),
                   paste(select,c("FC","pvalue","p.adj"),sep = "_"))
    sub_res = d[,select_col]
    sub_res[,select_col[5:8]] %<>% lapply(function(x) as.numeric(x))
    
    # filtering1
    first_filter = sub_res[,paste0(select,"_FC")] > 0 & sub_res[,paste0(select,"_pvalue")] < 0.05
    print(table(first_filter))
    connect = sub_res[first_filter,] %>% data.table::data.table()
    
    connect = connect[,c("R-L pair","Cell-cell pair")]
    connect = connect[,list(count=.N),by=`Cell-cell pair`]
    print(dim(connect))
    
    # filtering2
    second_filter = sub_res[,paste0(select,"_FC")] >= FC & sub_res[,paste0(select,"_p.adj")] < 0.05
    print(table(second_filter))
    highly_connected_pairs = unique(sub_res[second_filter,"Cell-cell pair"])
    print(length(highly_connected_pairs))

    
    # prepare cell-cell pair with 0 count
    All_cell_pair = paste0(rep(df$cell.types, each = 59), "_",rep(df$cell.types, time = 59))
    
    if(count_0){
        table(All_cell_pair %in% connect$`Cell-cell pair`)
        connect1 = data.frame(`Cell-cell pair` = All_cell_pair[!(All_cell_pair %in% connect$`Cell-cell pair`)],
                          `count` = 0)
        colnames(connect1)[1] = "Cell-cell pair"
        connect %<>% as.data.frame %>% rbind(connect1)
    }
    
    connect %<>% as.data.frame %>% cbind(str_split(string = connect$`Cell-cell pair`,pattern = "_",simplify = T))
    table(keep <- (connect$`1` %in% select_celltypes) & (connect$`2` %in% select_celltypes))
    connect = connect[keep,]
    #connect = connect[1:300,]
    print(dim(connect))
    connect$`1` %<>% gdata::reorder.factor(new.order=select_celltypes)
    connect %<>%   arrange(`1`)
    connect$value = connect$count/max(connect$count)
    ##==========  Hierarchical edge bundling ==================
    # Transform the adjacency matrix in a long format
    colnames(connect) = c("pair","count","from","to","value")
    connect = connect[,c("from","to","count","value","pair")]
    # Number of connection per person
    as_tibble(connect) %>%
        group_by(to) %>%
        summarize(sum=sum(count)) -> vertices
    colnames(vertices) <- c("name", "sum")
    vertices$value = vertices$sum/max(vertices$sum)
    vertices = vertices[,c("name", "value","sum")]

    # assign group
    group <- df[match(vertices$name, df$cell.types),"Cell.group"]
    group = na.omit(group$Cell.group)
    #Reorder dataset and make the graph
    
    vertices <- vertices %>% 
        mutate( group = group) %>%
        arrange(factor(name, levels = cell_order))

    # Add label angle
    number_of_bar=nrow(vertices)
    vertices$id = seq(1, nrow(vertices))
    angle= 360 * (vertices$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    vertices$hjust <- ifelse(angle > 90 & angle<270, 1, 0)
    vertices$angle <- ifelse(angle > 90 & angle<270, angle+180, angle)
    
    
    #dim(connect)
    # The connection object must refer to the ids of the leaves:
    edges = connect[connect$pair %in% highly_connected_pairs,c("from","to")]
    # remove link to self
    Remove = as.character(edges$from) == as.character(edges$to)
    if(any(Remove)) edges = edges[!Remove,]
    rownames(edges) = 1:nrow(edges)
    from  <-  match(as.character(edges$from), as.character(vertices$name))
    to  <-  match(as.character(edges$to), as.character(vertices$name))
    
    
    # prepare a vector of n color in the viridis scale
    mycolor <- unique(df$hexcolor)[match(unique(df$Cell.group),unique(vertices$group))]
    # Create a graph object with igraph
    mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )
    
    # Make the graph
ggraph(mygraph, layout = 'circle') +  
    #geom_edge_link(edge_colour="black",edge_alpha=1,edge_width=1)+
    geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, width=0.9, aes(colour=group), tension = 0.9) +
    #scale_edge_colour_distiller(palette = "RdPu") +
        geom_node_text(aes(label=paste("  ",name,"  "), angle=angle, hjust=hjust, colour=group), size=7) +
        geom_node_point(aes(size=value, color=group, fill=group,alpha=0.5)) +
        scale_size_continuous( range = c(0.1,10) ) +
        scale_color_manual(values=mycolor) +
        theme_void() +
        theme(
            legend.position="none",
            plot.margin=unit(c(0,0,0,0), "cm"),
            panel.spacing=unit(c(0,0,0,0), "cm")
        ) +
        expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))
    
#    jpeg(paste0(save.path,"new_chordDiagram_",select,".jpeg"), family = "Arial",
#         units="in", width=7, height=7,res=900)
#    print(g)
#    dev.off()
    }
Hierarchical_Edge_Bundling(d = df_res, select ="P vs D+T", FC = expm1(2), count_0 = T,
                           save.path = path)








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
