setwd('/Users/weizhao/Documents/NeuronChat-ST/MERFISH/processed_data/')
# locations and cell type annotations of MERFISH for MOp 
loc_centroid_subset <- readRDS('meta_MERFISH.rds')
# scRNA-seq data for MOp, with cell type annotations
MOp_sc_list <- readRDS(file = 'MOp_sc_list.rds')
MOp_sc_list$meta_sc$subclass_label <- sub(' CTX','',MOp_sc_list$meta_sc$subclass_label)
cell_type_use <- unique(MOp_sc_list$meta_sc$subclass_label)
cell_type_use_merfish <- unique(loc_centroid_subset$subclass)

# cell types shared by MERFISH & scRNA-seq
cell_type_shared <- sort(intersect(cell_type_use_merfish,cell_type_use),method='radix')
cell_type_shared <- cell_type_shared[8:12] ## i.e., c("Lamp5","Pvalb","Sncg","Sst","Vip")

# subsetting MERFISH meta data for single slice ('mouse1_slice212') and cell_type_shared
meta_single_slice <- loc_centroid_subset[loc_centroid_subset$slice_id=='mouse1_slice212' & loc_centroid_subset$subclass %in% cell_type_shared,]
# calculating cell proximity enrichment score for single slice 
cell_proximity <- cell_proximity_enrichment_score_single(meta_single_slice,celltype_label='subclass',centroid_x='centroid_x', centroid_y='centroid_y', 
                                       thresh_dist=400,permutation_number=1000)

# subsetting MERFISH meta data for cell_type_shared but all slices
meta_all_slices  <- loc_centroid_subset[loc_centroid_subset$subclass %in% cell_type_shared,]
# calculating cell proximity enrichment score for multiple slices 
cell_proximity <- cell_proximity_enrichment_score_multiple(meta_all_slices,celltype_label='subclass',centroid_x='centroid_x', centroid_y='centroid_y', slice_id='slice_id',
                                       thresh_dist=400,permutation_number=1000)

# circle plot of cell proximity 
library(scales);colours <- hue_pal()(length(cell_type_shared)+1);  colours[3] <- 'purple';colours[4] <- 'green'; 
netVisual_proximity(cell_proximity$log2CPscore_mtx,color.use = colours[1:5])
# barplot of cell proximity
barplot_proximity(cell_proximity$CPScore_df)

# NeuronChat analysis
GABA_idx <- which(MOp_sc_list$meta_sc$subclass_label %in% cell_type_shared);
MOp_sc_GABA <- MOp_sc_list$MOp_sc[GABA_idx,]
meta_sc_GABA <- MOp_sc_list$meta_sc[GABA_idx,]
x <- createNeuronChat(t(MOp_sc_GABA),DB='mouse',group.by = meta_sc_GABA$subclass_label);
x <- run_NeuronChat(x,M=100)
net_aggregated_x <- net_aggregation(x@net,method = 'weight')
# original network
netVisual_circle_neuron(net_aggregated_x, color.use = colours[1:5],vertex.size.max = 10,vertex.label.cex=1.5,edge.width.max=10)
# spatially constrained network (by filtering out links with cell proximity enrichment scores (stored in cell_proximity$log2CPscore_mtx as a matrix form) lower than 0)
netVisual_circle_neuron(net_aggregated_x*(cell_proximity$log2CPscore_mtx>0), color.use = colours[1:5],vertex.size.max = 10,vertex.label.cex=1.5,edge.width.max=10)

# spatialMap & spatial distribution 
p1 <- spatialMap(meta_single_slice,celltype_label='subclass',centroid_x='centroid_x', centroid_y='centroid_y', pt.size=5,pt.alpha=0.8,font.size=20,legend.symbol.size=5)
p2 <- spatialDistribution(meta_single_slice, celltype_label='subclass',centroid='centroid_y', curve.alpha=0.5,font.size=20,legend.symbol.size=5)
p1+p2
