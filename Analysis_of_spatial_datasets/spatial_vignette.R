setwd("/Users/weizhao/Documents/NeuronChat-ST/ST_datasets_processed/")
spaMap <- function(meta,point.size=1,alpha=0.8){
  title_x <- c('spatial map');attr_x <- c('cell_type')
  p2 <- ggplot2::ggplot() + theme_bw()
  p2 <- p2+ ggplot2::geom_point(data = meta[order(meta[,attr_x],decreasing = F),], aes_string(x = 'sdimx', y = 'sdimy',colour = attr_x),
                                show.legend = T, size = point.size, alpha =alpha) + scale_color_manual(values=colours) + labs(title_x)+
    guides(color = guide_legend(override.aes = list(size=5))) +  
    theme(text=element_text(size=20), #change font size of all text
          axis.text=element_text(size=20), #change font size of axis text
          axis.title=element_text(size=20), #change font size of axis titles
          plot.title=element_text(size=20,hjust = 0.5), #change font size of plot title
          legend.text=element_text(size=20), #change font size of legend text
          legend.title=element_text(size=20)) #change font size of legend title 
  return(p2)
}
## MERFISH, 2021, https://www.nature.com/articles/s41586-021-03705-x
load("/Users/weizhao/Documents/NeuronChat-ST/ST_datasets_processed/MERFISH/visual_spatial_MERFISH.rda")
# set.seed(123)
# x <- createNeuronChat(t(mat_sc),DB='mouse',group.by = as.character(meta_sc$subclass_label)) # input data should be gene by cell matrix
# x <- run_NeuronChat(x,M=100)
# net_aggregated_x <- net_aggregation(x@net,method = 'weight')
# net <- net_aggregated_x/max(net_aggregated_x)
# # only keep top10 links with highest total weight
# net_sort <- sort(net,index.return	=T,decreasing = T)
# net_top10 <-net*0;net_top10[net_sort$ix[1:10]] <-net[net_sort$ix[1:10]]
# save(mat_sc,meta_sc, meta, net, net_top10,file='visual_spatial_MERFISH.rda')

# circle plot of communication network
library(scales);colours <- hue_pal()(dim(net)[1]);
netVisual_circle_neuron(net,edge.width.max = 10,color.use = colours,vertex.label.cex = 1.5)
netVisual_circle_neuron(net_top10,edge.width.max = 10,color.use = colours,vertex.label.cex = 1.5)
# spatial map of cell types
spaMap(meta) + labs(x ="x coordinate", y = "y coordinate",color='Glutamatergic \n subclass')

## seqFISH, 2019, https://www.nature.com/articles/s41586-019-1049-y
load("/Users/weizhao/Documents/NeuronChat-ST/ST_datasets_processed/SeqFISH/visual_spatial_seqFISH.rda")
# set.seed(123)
# x <- createNeuronChat(mat,DB='mouse',group.by = as.character(meta$cell_type))
# x <- run_NeuronChat(x,M=100,strict = 0) # use strict=0 to allow calculation for interaction pairs with missing genes, e.g., 'Gls' is missing from this dataset
# net_aggregated_x <- net_aggregation(x@net,method = 'weight')
# net <- net_aggregated_x/max(net_aggregated_x)
# # only keep top10 links with highest total weight
# net_sort <- sort(net,index.return	=T,decreasing = T)
# net_top10 <-net*0;net_top10[net_sort$ix[1:10]] <-net[net_sort$ix[1:10]]
# save(mat,meta, net, net_top10,file='visual_spatial_seqFISH.rda')

# circle plot of communication network
library(scales);colours <- hue_pal()(dim(net)[1]);
netVisual_circle_neuron(net,edge.width.max = 10,color.use = colours,vertex.label.cex = 1.5)
netVisual_circle_neuron(net_top10,edge.width.max = 10,color.use = colours,vertex.label.cex = 1.5)
# spatial map of cell types
spaMap(meta,point.size = 5) + labs(x ="x coordinate", y = "y coordinate",color='cell type')


## Visium, 2020, https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain
load("/Users/weizhao/Documents/NeuronChat-ST/ST_datasets_processed/Visium/visual_spatial_Visium.rda")
set.seed(123)
# x <- createNeuronChat(mat,DB='mouse',group.by = as.character(meta$cell_type))
# x <- run_NeuronChat(x,M=100)
# net_aggregated_x <- net_aggregation(x@net,method = 'weight')
# net <- net_aggregated_x/max(net_aggregated_x)
# net <- net[levels(meta$cell_type),levels(meta$cell_type)]
# net_sort <- sort(net,index.return	=T,decreasing = T)
# net_top10 <-net*0;net_top10[net_sort$ix[1:10]] <-net[net_sort$ix[1:10]]
library(scales);colours <- hue_pal()(dim(net)[1]);colours <- colours[c(3,2,1,4,5,6,7)] 
netVisual_circle_neuron(net,edge.width.max = 10,color.use = colours,vertex.label.cex = 1.5)
netVisual_circle_neuron(net_top10,edge.width.max = 10,color.use = colours,vertex.label.cex = 1.5)
# save(mat,meta, net, net_top10,file='visual_spatial_Visium.rda')
spaMap(meta,point.size=2,alpha = 1) + labs(x ="x coordinate", y = "y coordinate",color='spot cluster')