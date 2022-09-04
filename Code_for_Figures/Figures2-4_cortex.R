library(NeuronChat)
library(data.table)
library(CellChat)
library(patchwork)
library(rhdf5)
library(bigstatsr)
library(Matrix)
library(pryr);library(dplyr);
library(gtools)
library(visNetwork)
library(RColorBrewer)
library(ComplexHeatmap)
library(PRROC)
options(stringsAsFactors = FALSE)

###### data preparation
#######################
### target_df for VISp and ALM
setwd('/Users/weizhao/Documents/CellChat/Rcode_cellchat_VISp_ALM/')
rna <- fread('GSE115746_cells_exon_counts.csv.gz') # # data from Tasic 2018 Nature
rna_new <- data.frame(rna[,2:ncol(rna)],row.names=rna$V1)
meta_data <- fread('GSE115746_complete_metadata_28706-cells.csv')
# <- meta_data[meta_data$source_name=='Primary Visual Cortex (VISp)',]
meta_use <- meta_data[meta_data$source_name %in% c("Anterior Lateral Motor Cortex (ALM)","Primary Visual Cortex (VISp)"),]
meta_use <- meta_use[meta_use$cell_class %in% c("Glutamatergic"),]#"Glutamatergic","Non-Neuronal","Endothelial"),]
meta_use <- meta_use[meta_use$cell_subclass != "",]
meta_use <- meta_use[meta_use$cell_subclass != 'CR',]
cell_use <- meta_use$sample_name
mat <- rna_new[,cell_use[cell_use %in% names(rna_new)]]
meta <- meta_use[meta_use$sample_name %in% names(rna_new),]
row.names(meta) <- meta$sample_name

# normalization
library.size <- Matrix::colSums(mat)
mat_tmp <- t(t(mat)/ library.size) * 10000
mat_norm <- log1p(mat_tmp)
target_mat <- mat_norm[target_genes[which(target_genes %in% rownames(mat_norm))],]
target_df <- data.frame(t(target_mat))
target_df$cell_subclass <- meta$cell_subclass # (meta$cluster_label)
VISp_idx <- which(meta$source_name=='Primary Visual Cortex (VISp)')
ALM_idx <- which(meta$source_name=='Anterior Lateral Motor Cortex (ALM)')
target_df_VISp <- target_df[VISp_idx,]
target_df_ALM <- target_df[ALM_idx,]
# saveRDS(target_df_ALM,'target_df_ALM.rds')
# saveRDS(target_df_VISp,'target_df_VISp.rds')


setwd('/Users/weizhao/Documents/CellChat/Rcode_cellchat_cortex_hippocampus/CTX-HPF_deposit/')
filename_data <- 'expression_matrix_SSv4.hdf5' # data from Yao 2021 Cell
h5ls(filename_data)
i <- as.integer(h5read(filename_data,"/data/exon/i"))
p <- as.integer(h5read(filename_data,"/data/exon/p"));dp <- diff(p)
x <- as.integer(h5read(filename_data,"/data/exon/x"))
rna <- sparseMatrix(i=i, p=p,x=x,dims = c(73363,45768),index1=FALSE)
gene_symbol <- h5read(filename_data,"/gene_names")
sample_names <- h5read(filename_data,"/sample_names")
dimnames(rna) <- list(sample_names,gene_symbol)
meta_data <- fread('metadata_ssv4.csv')
meta_use <- meta_data[ which(meta_data$sample_name %in% sample_names),]

regions_all <- names(table(meta_use$region_label))
cortex_regions <-  c("ACA","AI","ALM","AUD","CLA","GU","MOp","ORB","PL-ILA","PTLp","RSP","RSPv","SSp","SSs","TEa-PERI-ECT","VIS","VISp")
HPF_regions <- c("ENTl","ENTm","HIP","PAR-POST-PRE","SUB-ProS")

cell_use <- meta_use$sample_name
mat <- rna[cell_use,]
meta <- meta_use;rm(i);rm(p);rm(x);rm(rna);gc()
# aggregate((mat[,'Gad2']),list(meta$subclass_label),FUN=triMean)

# normalization
library.size <- Matrix::rowSums(mat)
mat_tmp <- (mat/ library.size) * 10000
mat_norm <- log1p(mat_tmp)
rm(mat_tmp);gc()
# saveRDS(mat_norm,file='mat_cortex.rds')
setwd('/Users/weizhao/Documents/CellChat/Rcode_cellchat_VISp_ALM/')
interactionDB <- readRDS('interactionDB_mouse.rds')
target_genes <- c()
for(j in 1:length(interactionDB)){
  target_genes <- c(target_genes,interactionDB[[j]]$lig_contributor,interactionDB[[j]]$receptor_subunit,interactionDB[[j]]$targets_up,interactionDB[[j]]$targets_down,interactionDB[[j]]$activator,interactionDB[[j]]$inhibitor, interactionDB[[j]]$interactors)
}
target_genes <- unique(target_genes)
target_mat <- mat_norm[,target_genes[which(target_genes %in% colnames(mat_norm))]]
target_df <- as.data.frame(as.matrix(target_mat))

#### target_df for all cortex regions
target_df$cell_subclass <- meta$subclass_label
# saveRDS(target_df,file='target_df_cortex.rds')

### target_df for VISp projection target regions
VISp_proj_idx <- which( (meta$region_label %in% c('VIS','VISp','ACA','RSP','AUD')|  (meta$region_label %in% c('HIP') & meta$subclass_label=='DG')) & meta$class_label=='Glutamatergic')
table(meta$region_label[VISp_proj_idx])
target_df_VISp_proj <- target_df[VISp_proj_idx,]
target_df_VISp_proj$cell_subclass <- meta$region_label[VISp_proj_idx]
target_df_VISp_proj$cell_subclass[target_df_VISp_proj$cell_subclass=='VISp'] <- 'VISp-c'
target_df_VISp_proj$cell_subclass[target_df_VISp_proj$cell_subclass=='HIP'] <- 'HIP_DG'
target_df_VISp_proj$cell_subclass[target_df_VISp_proj$cell_subclass!='HIP_DG'] <- paste('CTX_',target_df_VISp_proj$cell_subclass[target_df_VISp_proj$cell_subclass!='HIP'],sep = '')
# saveRDS(target_df_VISp_proj,file='target_df_VISp_proj.rds')
target_df_VISp_proj$cell_subclass[target_df_VISp_proj$cell_subclass=='CTX_HIP_DG'] <- 'HIP_DG'

### target_df for ALM projection target regions
ALM_proj_idx <- which(meta$region_label %in% c('MOp','ALM','ORB','RSP','SSp','SSs') & meta$class_label=='Glutamatergic')
table(meta$region_label[ALM_proj_idx])
target_df_ALM_proj <- target_df[ALM_proj_idx,]
target_df_ALM_proj$cell_subclass <- meta$region_label[ALM_proj_idx]
target_df_ALM_proj$cell_subclass[target_df_ALM_proj$cell_subclass=='ALM'] <- 'ALM-c'
target_df_ALM_proj$cell_subclass[target_df_ALM_proj$cell_subclass=='ORB'] <- 'ORB-c'
target_df_ALM_proj$cell_subclass[target_df_ALM_proj$cell_subclass!='HIP'] <- paste('CTX_',target_df_ALM_proj$cell_subclass[target_df_ALM_proj$cell_subclass!='HIP'],sep = '')
# saveRDS(target_df_ALM_proj,file='target_df_ALM_proj.rds')


setwd('/Users/weizhao/Documents/CellChat/Rcode_cellchat_VISp_ALM/')
target_df <- readRDS(file='target_df_cortex.rds')
choose_single <- function(region_name,cell_class=NULL){
  if(is.null(cell_class)){cell_class <- names(table(meta$class_label))}
  cell_idx <- which(meta$region_label %in% region_name & meta$class_label %in% cell_class & meta$subclass_label %in% subclass_CTX)
  target_df_tmp <- target_df[cell_idx,]
  return(target_df_tmp)
}
cortex_cellnumber <- table(meta$region_label)[cortex_regions]
cortex_1K <- names(cortex_cellnumber[cortex_cellnumber>1000])
cortex_1K <- c('ALM','VISp')
subclass_CTX <- names(table(meta$subclass_label[meta$class_label=='Glutamatergic']))[grep('CTX',names(table(meta$subclass_label[meta$class_label=='Glutamatergic'])))]
set.seed(1234)
cortex_list <- lapply(cortex_1K, function(y){
  #target_df_single <- choose_single(y,'Glutamatergic')
  region_name <- y; cell_class <- 'Glutamatergic';if(is.null(cell_class)){cell_class <- names(table(meta$class_label))}
  cell_idx <- which(meta$region_label %in% region_name & meta$class_label %in% cell_class & meta$subclass_label %in% subclass_CTX)
  target_df_single  <- target_df[cell_idx,]
  meta_tmp <- meta[cell_idx,];rownames(meta_tmp) <- meta_tmp$sample_name
  x <- createNeuronChat( t(as.matrix(target_df_single[,1:(dim(target_df_single)[2]-1)])),DB='mouse',group.by = target_df_single$cell_subclass,meta=meta_tmp);
  x <- run_NeuronChat(x,M=100)
  return(x)
})
names(cortex_list) <- cortex_1K

#### Figure 4a, 4b
par(mfrow=c(1,2))
for(j in c(1,2)){
  netVisual_circle_neuron(cortex_list[[j]]@net$Glu_Grin3a, title.name = paste('Glu_Grin3a -- ',names(cortex_list)[j]),arrow.size = 0.5,margin=0.3,edge.width.max=8)
}
par(mfrow=c(1,2))
for(j in c(1,2)){
  net_aggregated_x <- net_aggregation(cortex_list[[j]]@net,method='weight')
  netVisual_circle_neuron(net_aggregated_x, title.name = names(cortex_list)[j],arrow.size = 0.5,margin=0.3,edge.width.max=8)
}

#### Figure 4d
neuronchat_list <- mergeNeuronChat(cortex_list, add.names = names(cortex_list))
neuronchat_list <- computeNetSimilarityPairwise_Neuron(neuronchat_list, slot.name = "net", type = "functional",comparison = c(1,2))
neuronchat_list <- netEmbedding(neuronchat_list,slot.name = "net_analysis", type = "functional",comparison = c(1,2))
neuronchat_list <- netClustering(neuronchat_list, slot.name = "net_analysis", type = "functional",comparison = c(1,2),k = 5)
netVisual_embeddingPairwise_Neuron(neuronchat_list, slot.name = "net_analysis", type = "functional", label.size = 3.5,comparison=c(1,2),pathway.remove.show = FALSE,pathway.labeled = F)
netVisual_embeddingPairwiseZoomIn_Neuron(neuronchat_list, slot.name = "net_analysis", type = "functional", label.size = 5,comparison=c(1,2),nCol=2)

#### Figure 4c
g1 <- rankNet_Neuron(neuronchat_list,mode='comparison',measure = c("count"),comparison = 1:2,do.stat = F,tol = 0.1,stacked = F,font.size = 11)
g2 <- rankNet_Neuron(neuronchat_list,mode='comparison',measure = c("weight"),comparison = 1:2,do.stat = F,tol = 0.1,stacked = F,font.size = 11)#+scale_y_continuous(trans = "log1p",breaks = c(0,1,10,50,200),labels =c(0,1,10,50,200)) #scale_y_continuous(breaks = log2(c(1,10,120)), labels =log2(c(1,10,120))) ## + scale_y_log10()
g1+g2
# compareInteractions_Neuron(neuronchat_list,measure = c("weight"),comparison = c(3,15),group=c(1,2),show.legend = F)
# compareInteractions_Neuron(neuronchat_list,measure = c("count"),comparison = c(3,15),group=c(1,2),show.legend = F )
cluster_info <- neuronchat_list@net_analysis$similarity$functional$group$`1-2`
cluster_info_region <- cluster_info;names(cluster_info_region) <- sub(".*--", "", names(cluster_info_region))
cluster_info_pathway <- cluster_info;names(cluster_info_pathway) <- sub("--.*","", names(cluster_info_pathway))
df_info <- data.frame(cluster=cluster_info,region=names(cluster_info_region),pathway=names(cluster_info_pathway))
df_info <- df_info[order(df_info$pathway),]

#### Figure 4e
net3_15 <- neuronchat_list@net[c(3,15)]
net3 <- net3_15[[1]];names(net3) <- paste(names(net3),'--ALM',sep='')
net15 <- net3_15[[2]];names(net15) <- paste(names(net15),'--VISp',sep='')
net3_15_list <- append(net3,net15)
interaction_group <- neuronchat_list@net_analysis$similarity$functional$group$`3-15`
hlist <- list();gb_heatmap <- list()
grid.newpage(); x_seq <- c(0,0.2,0.4,0.6,0.8)
for(j in 1:length(unique(interaction_group))){
  net_aggregated_group2 <- net_aggregation(net3_15_list[names(interaction_group[interaction_group==j])],method = 'weight')
  library(RColorBrewer);col_map = brewer.pal(8,"YlOrBr");
  h <- Heatmap(net_aggregated_group2, name = "Weight",
                        col = col_map,
                        cluster_rows = FALSE,cluster_columns=FALSE,
                        row_names_side='left',column_names_side='bottom',
                        #row_title='Sender',row_title_side='left',
               row_title_gp = gpar(fontsize = 16),
                        column_title='Receiver',column_title_side = "bottom",column_title_gp = gpar(fontsize = 16),column_names_rot = 60)
  gb_heatmap[[j]] = grid.grabExpr(draw(h, padding = unit(c(2, 2, 2, 2), "mm")) )
  pushViewport(viewport(x = x_seq[j], y = 1, width = 0.19, height = 1, just = c("left", "top"),xscale = c(0, 1), yscale = c(0, 1)));grid.draw(gb_heatmap[[j]]);popViewport()
}


#####################
#####################
####### VISp, Figure 3
y <- 'VISp'
VISp <- lapply('VISp',function(y){
  region_name <- y; cell_class=NULL;if(is.null(cell_class)){cell_class <- names(table(meta$class_label))}
  cell_idx <- which(meta$region_label %in% region_name & meta$class_label %in% cell_class & !(meta$subclass_label %in%c('Car3','CR','DG','L2/3 IT PPP','L5/6 IT TPE-ENT')))
  target_df_single  <- target_df[cell_idx,]
  meta_tmp <- meta[cell_idx,];rownames(meta_tmp) <- meta_tmp$sample_name
  df_group <- meta_tmp[!duplicated(meta_tmp$subclass_label),c('class_label','subclass_label')]
  group <- structure(df_group$class_label,names=df_group$subclass_label)
  #group <- group[!(names(group) %in% c('Car3','CR','DG','L2/3 IT PPP','L5/6 IT TPE-ENT'))]
  x <- createNeuronChat( t(as.matrix(target_df_single[,1:(dim(target_df_single)[2]-1)])),DB='mouse',group.by = target_df_single$cell_subclass,meta=meta_tmp);
  x <- run_NeuronChat(x,M=100)
  net_aggregated_x <- net_aggregation(x@net,method = 'weight')
  # Figure3a, circle plot
  netVisual_circle_neuron(net_aggregated_x, title.name = y,group=group, vertex.weight = rowSums(net_aggregated_x))
  #netVisual_circle_zw(net_aggregated_x, title.name = y,group=group, vertex.weight = rowSums(net_aggregated_x))

  # # Figure3a, chordDiagram
  CellChat::netVisual_chord_cell_internal(net_aggregated_x, group = group,lab.cex=1.3)
  ## Figure 3a Heatmap
  library(RColorBrewer);par(mfrow=c(1,1));col_map = brewer.pal(8,"YlOrBr"); ng <- length(unique(group));
  left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill =scPalette(ng)),labels = sort(unique(group)),labels_gp = gpar(col = "white")))
  bottom_Annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill =scPalette(ng)), labels = sort(unique(group)),labels_gp = gpar(col = "white")))
  h1 <- ComplexHeatmap::Heatmap(net_aggregated_x, name = "Weight", left_annotation =left_Annotation,bottom_annotation = bottom_Annotation,
                #top_annotation = column_ha, right_annotation = row_ha,
                col = col_map,
                cluster_rows = FALSE,cluster_columns=FALSE,
                split=as.character(group[rownames(net_aggregated_x)]),column_split = as.character(group[rownames(net_aggregated_x)]),
                row_names_side='left',column_names_side='bottom',
                row_title='Sender',row_title_side='left',row_title_gp = gpar(fontsize = 16),
                column_title='Receiver',column_title_side = "bottom",column_title_gp = gpar(fontsize = 16),column_names_rot = 60
  )#,
  draw(h1)
  heatmap_single(x,interaction_name='Glu_Gria2',group=group)
  return(x)
})

### Figure 3b
VISp <- list()
VISp[[1]] <- x
selectK_Neuron(VISp[[1]],pattern = "outgoing")
selectK_Neuron(VISp[[1]],pattern = "incoming")
VISp[[1]]<- identifyCommunicationPatterns_Neuron(VISp[[1]], slot.name = "net", pattern = c("outgoing"), k=4,height = 18, thresh_quantile = 0)
VISp[[1]] <- identifyCommunicationPatterns_Neuron(VISp[[1]], slot.name = "net", pattern = c("incoming"), k=4,height = 18)
gg1 <- netAnalysis_river_Neuron(VISp[[1]],slot.name = "net", pattern = c("outgoing"),font.size = 2.5,cutoff.1 = 0.5,cutoff.2=0.5,top.n1=1e4,top.n2=1e4)
gg2 <- netAnalysis_river_Neuron(VISp[[1]],slot.name = "net", pattern = c("incoming"),font.size = 2.5,cutoff.1 = 0.5,cutoff.2=0.5,top.n1=1e4,top.n2=1e4)
gg1
gg2

library(ggalluvial)
VISp_tmp <- VISp[[1]];names(VISp_tmp@info) <- names(VISp_tmp@net)
top50 <- names(sort(VISp_tmp@info,decreasing = T)[1:60]);top50 <- names(VISp_tmp@info)[names(VISp_tmp@info) %in% top50]
VISp_tmp@net <- VISp_tmp@net[top50]
VISp_tmp <- identifyCommunicationPatterns_Neuron(VISp_tmp, slot.name = "net", pattern = c("outgoing"), k=4,height = 18, thresh_quantile = 0)
VISp_tmp <- identifyCommunicationPatterns_Neuron(VISp_tmp, slot.name = "net", pattern = c("incoming"), k=4,height = 18)
gg1 <- netAnalysis_river_Neuron(VISp_tmp,slot.name = "net", pattern = c("outgoing"),font.size = 3,cutoff.1 = 0.5,cutoff.2=0.5,top.n1=1e4,top.n2=1e4)
gg2 <- netAnalysis_river_Neuron(VISp_tmp,slot.name = "net", pattern = c("incoming"),font.size = 3,cutoff.1 = 0.5,cutoff.2=0.5,top.n1=1e4,top.n2=1e4)
gg1
gg2


### supplementary figure related to Figure 3
g1 <- rankNet_Neuron(VISp[[1]],slot.name = "net",measure = c("weight"),mode='single',font.size = 5)
g2 <- rankNet_Neuron(VISp[[1]],slot.name = "net",measure = c("count"),mode='single',font.size = 5)
g1+g2


### Figures 3c & 3d
VISp[[1]] <- computeNetSimilarity_Neuron(VISp[[1]],type='functional')
VISp[[1]]  <- netEmbedding(VISp[[1]], slot.name = "net_analysis", type = "functional")
VISp[[1]] <- netClustering(VISp[[1]],type='functional',slot.name = "net_analysis",k=5)
netVisual_embedding_Neuron(VISp[[1]], type = "functional", label.size = 5,pathway.remove.show = F)
netVisual_embeddingZoomIn_Neuron(VISp[[1]], type = "functional", nCol = 2,label.size = 3)

### supplementary figure related to Figure 3
hlist <- list();gb_heatmap <- list()
grid.newpage(); x_seq <- c(0, 1/3,2/3);y_sed <- c(1,0.5)
position_df <- as.matrix(expand.grid(x_seq,y_sed))
for(j in 1:5){
net_aggregated_group2 <- net_aggregation(x@net[names(VISp[[1]]@net_analysis$similarity$functional$group$single[VISp[[1]]@net_analysis$similarity$functional$group$single==j])],method = 'weight')
library(RColorBrewer);par(mfrow=c(1,1));col_map = brewer.pal(8,"YlOrBr"); ng <- length(unique(group));
left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill =scPalette(ng)),labels = sort(unique(group)),labels_gp = gpar(col = "white")))
bottom_Annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill =scPalette(ng)), labels = sort(unique(group)),labels_gp = gpar(col = "white")))
h <- Heatmap(net_aggregated_group2, name = "Weight", left_annotation =left_Annotation,bottom_annotation = bottom_Annotation,
              #top_annotation = column_ha, right_annotation = row_ha,
              col = col_map,
              cluster_rows = FALSE,cluster_columns=FALSE,
              split=as.character(group[rownames(net_aggregated_group2)]),column_split = as.character(group[rownames(net_aggregated_group2)]),
              row_names_side='left',column_names_side='bottom',
              row_title='Sender',row_title_side='left',row_title_gp = gpar(fontsize = 16),
              column_title='Receiver',column_title_side = "bottom",column_title_gp = gpar(fontsize = 16),column_names_rot = 60
)#,
gb_heatmap[[j]] = grid.grabExpr(draw(h, column_title=paste('Group',j,sep = ' '),padding = unit(c(2, 2, 2, 2), "mm")) )
pushViewport(viewport(x = position_df[j,1], y = position_df[j,2], width = 0.32, height = 0.48, just = c("left", "top"),xscale = c(0, 1), yscale = c(0, 1)));grid.draw(gb_heatmap[[j]]);popViewport()
}

##### Figure2  for VISp projections
setwd('/Users/weizhao/Documents/CellChat/Rcode_cellchat_VISp_ALM/')
source("~/Documents/CellChat/Rcode_cellchat_VISp_ALM/0220/Projections (1).r", encoding = 'UTF-8')
sender.names <- sort(sender.names,method='radix')
receiver.names <- sort(receiver.names,method='radix')
### load two datasets
set.seed(123) ## for VISp then ALM
# for(j in 1:10){
# set.seed(123456)
# set.seed(j)
# df_sender <- readRDS('target_df_VISp.rds')
df_sender <- target_df_VISp
df_sender <- df_sender[df_sender$cell_subclass %in% sender.names,]
# df_receiver <- readRDS('target_df_VISp_proj.rds')
df_receiver <- target_df_VISp_proj
df_receiver <- df_receiver[df_receiver$cell_subclass %in% receiver.names,]
target_df <- rbind(df_sender,df_receiver)
tmp <- t(as.matrix(target_df[,1:(dim(target_df)[2]-1)]))
VISp <- createNeuronChat(tmp,DB='mouse',group.by = target_df$cell_subclass);
VISp <- run_NeuronChat(VISp,sender=sender.names, receiver=receiver.names,M=100)
net_aggregated_VISp <- net_aggregation(VISp@net,method='weight_threshold',cut_off = 0.8) # weight_threshold
net.names <- rownames(VISp_projection_mtx)
net_predicted <- VISp_projection_mtx;net_predicted[sender.names,receiver.names] <- net_aggregated_VISp[sender.names,receiver.names]
## Figure 2a, 2b
par(mfrow=c(1,2))
netVisual_circle_neuron(VISp_projection_mtx,group=group,vertex.label.cex=1.5,arrow.size = 0.8,margin = 0.4)
netVisual_circle_compare(VISp_projection_mtx,(net_predicted[net.names,net.names]/max(net_predicted[net.names,net.names])>0.028)*1,group=group,vertex.label.cex=1.5,arrow.size = 0.8,margin = 0.4)
## Figure 2c, 2d
library(PRROC)
PRROC_obj <- roc.curve(scores.class0 = c(net_predicted[sender.names,receiver.names]/max(net_predicted[net.names,net.names])), weights.class0=c(VISp_projection_mtx[sender.names,receiver.names]),curve=TRUE)
par(mar = c(5, 5, 4, 1),cex.lab=1.5, cex.axis=1.5,cex.main=1.5); plot(PRROC_obj,xlab = 'False positive rate')
PRROC_obj <- pr.curve(scores.class0 = c(net_predicted[sender.names,receiver.names]/max(net_predicted[net.names,net.names])), weights.class0=c(VISp_projection_mtx[sender.names,receiver.names]),curve=TRUE)
par(mar = c(5, 5, 4, 1),cex.lab=1.5, cex.axis=1.5,cex.main=1.5);plot(PRROC_obj)
### used for determine the threshold
accuracy<-c();sensitivity<-c();FPR<-c()
for(threshold in seq(0,1,1e-3)){
  predicted_net <- (net_predicted/max(net_predicted)>threshold)*1
  accuracy <- c(accuracy,1-sum(abs(predicted_net-VISp_projection_mtx))/length(receiver.names)/length(sender.names))
  sensitivity <- c(sensitivity,sum(c(predicted_net)*(c(VISp_projection_mtx)>0))/sum(c(VISp_projection_mtx)>0))
  FPR <- c(FPR,sum((c(predicted_net[sender.names,receiver.names])*(c(VISp_projection_mtx[sender.names,receiver.names])==0)))/sum((c(VISp_projection_mtx[sender.names,receiver.names])==0)))
}
print(max(accuracy))
seq(0,1,1e-3)[which(accuracy==max(accuracy))]

## Figure 2e, supplementary figure
## test AUROC for each interaction
VISp_label <- c(VISp_projection_mtx[sender.names,receiver.names])
auroc_all <- rep(0,length(VISp@net))
prroc_all <- rep(0,length(VISp@net))
for(j in 1:length(VISp@net)){
  if(max(VISp@net[[j]])==0){auroc_all[j]=0;prroc_all[j]=0} else{
    VISp_predict <- c(VISp@net[[j]][sender.names,receiver.names])/max(VISp@net[[j]])
    auroc_all[j] <- roc.curve(scores.class0 = VISp_predict, weights.class0=VISp_label,
                              curve=TRUE)$auc
    prroc_all[j] <- pr.curve(scores.class0 = VISp_predict, weights.class0=VISp_label,
                             curve=TRUE)$auc.integral
  }
}
names(auroc_all) <- names(VISp@DB)
auroc_all_sort <- sort(auroc_all,decreasing = T)
weight <- apply(simplify2array(VISp@net),3,sum)
count <- apply(simplify2array(VISp@net),3,FUN=function(x){sum(x>0)})
df_cor <- data.frame(AUROC=auroc_all,AUPRC=prroc_all,Info=weight,count=count);df_cor$interaction <- rownames(df_cor)
df_cor <- df_cor[df_cor$Info>0,]
t <- list()
cor.test(df_cor$AUROC,df_cor$AUPRC,method = 'spearman')
t[[1]] <- cor.test(df_cor$AUROC,df_cor$Info,method = 'spearman')
t[[2]] <- cor.test(df_cor$AUROC,df_cor$count,method = 'spearman')
t[[3]]<- cor.test(df_cor$AUPRC,df_cor$Info,method = 'spearman')
t[[4]] <- cor.test(df_cor$AUPRC,df_cor$count,method = 'spearman')
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
title_rho_p <- sapply(t,FUN=function(x){paste('rho = ',signif(x$estimate,2),'; p = ', signif(x$p.value, 2),sep = '')},simplify = T)
library(ggplot2)
library(hrbrthemes)
# mtcars dataset is natively available in R
# head(mtcars)
# A basic scatterplot with color depending on Species
p1 <- ggplot(df_cor, aes(x=AUROC, y=Info)) + labs(x ="AUROC", y = "Information flow",title=title_rho_p[1]) +geom_point(size=3)  + theme(text = element_text(size=20),axis.text.x = element_text(angle=0, hjust=1)) + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

p2 <- ggplot(df_cor, aes(x=AUROC, y=count)) + labs(title=title_rho_p[2],
                                                   x ="AUROC", y = "Count of links") +geom_point(size=3)  + theme(text = element_text(size=20),
                                                                                                                  axis.text.x = element_text(angle=0, hjust=1)) +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

p11 <- ggplot(df_cor, aes(x=AUPRC, y=Info)) + labs(title=title_rho_p[3],
                                                   x ="AUPRC", y = "Information flow") +geom_point(size=3)  + theme(text = element_text(size=20),
                                                                                                                    axis.text.x = element_text(angle=0, hjust=1)) +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

p22 <- ggplot(df_cor, aes(x=AUPRC, y=count)) + labs(title=title_rho_p[4],
                                                    x ="AUPRC", y = "Count of links") +geom_point(size=3)  + theme(text = element_text(size=20),
                                                                                                                   axis.text.x = element_text(angle=0, hjust=1)) +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

p1 + p11 +p2 + p22

p3 <- ggplot(df_cor,aes(x = Info,y = reorder(interaction,Info))) + geom_col(fill='orange') + theme(text = element_text(size=18)) + labs(y= "Interaction", x = "Information flow") +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
p4 <- ggplot(df_cor,aes(x = count,y = reorder(interaction,count))) + geom_col(fill='orange') + theme(text = element_text(size=18)) + labs(y= "Interaction", x = "Count of links") +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
p5 <- ggplot(df_cor,aes(x = AUROC,y = reorder(interaction,AUROC))) + geom_col(fill='orange') + theme(text = element_text(size=18)) + labs(y= "Interaction", x = "AUROC") +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) + geom_vline(xintercept=0.5, linetype="dashed", color = 'red', size=0.5)
p6 <- ggplot(df_cor,aes(x = AUPRC,y = reorder(interaction,AUPRC))) + geom_col(fill='orange') + theme(text = element_text(size=18)) + labs(y= "Interaction", x = "AUPRC") +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) + geom_vline(xintercept=sum(VISp_label)/length(sender.names)/length(receiver.names), linetype="dashed", color = 'red', size=0.5)
p5+p6+p3+p4+plot_layout(ncol = 4)

rownames(VISp@outgoing) <- rownames(VISp@net[[1]])



library(ggplot2)
df_auroc <- data.frame(interaction=names(auroc_all_sort),auroc=auroc_all_sort); df_auroc <- df_auroc[df_auroc$auroc!=0,]
ggplot(df_auroc,aes(x = auroc,y = reorder(interaction,auroc))) + geom_col(fill='orange') + theme(text = element_text(size=18)) + labs(y= "Interaction", x = "AUROC")+ geom_vline(xintercept=0.5, linetype="dashed", color = 'red', size=0.5)

lig_tar_name_plot <- rownames(df_cor)[order(df_cor$AUROC,decreasing = T)]
info_mat <- rbind(info=VISp@info,VISp@outgoing[sender.names,],VISp@incoming[receiver.names,])
column_ha = HeatmapAnnotation(Incoming = anno_barplot(VISp@incoming[receiver.names,interaction_name],border = F,height = unit(2,'cm'),gp = gpar(fill = "lightskyblue1")),annotation_name_side='left',annotation_name_rot = 0,annotation_label = 'Incoming strength')
row_ha = rowAnnotation(ligand = anno_barplot(VISp@outgoing[sender.names,interaction_name],border = F,width = unit(2,'cm'),gp = gpar(fill = "lightskyblue1")),annotation_label = 'ligand \n abundance')
Heatmap(info_mat[c(sender.names,receiver.names),lig_tar_name_plot], name = "Ligand \n Abundance",
        # top_annotation = column_ha, right_annotation = row_ha,
        cluster_rows = F,cluster_columns = F,column_names_rot = 45,
        row_names_side = 'left',col=col_map)
###
## barplot for ligand score & receptor score
col_map=colorspace::diverge_hsv(10)
ligand_profile <- VISp@ligand.abundance[sender.names,VISp@info>0]
ligand_profile <-  ligand_profile[,order(colSums(ligand_profile),decreasing = T)]
target_profile <- VISp@target.abundance[receiver.names,VISp@info>0]
target_profile <-  target_profile[,order(colSums(target_profile),decreasing = T)]
h1 <- Heatmap(ligand_profile,
              name = "Ligand \n Abundance", cluster_rows = F,cluster_columns = F,column_names_rot = 60,
              row_names_side = 'left',col=col_map,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 12))
h2 <- Heatmap(target_profile,
              name = "Target \n Abundance", cluster_rows = F,cluster_columns = F,column_names_rot = 60,
              row_names_side = 'left',col=col_map,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 12))

h1 <- Heatmap(t(ligand_profile),
              name = "Ligand \n Abundance", cluster_rows = F,cluster_columns = F,column_names_rot = 60,
              row_names_side = 'left',col=col_map,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 12))
h2 <- Heatmap(t(target_profile),
              name = "Target \n Abundance", cluster_rows = F,cluster_columns = F,column_names_rot = 60,
              row_names_side = 'left',col=col_map,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 12))
grid.newpage()
pushViewport(viewport(x = 0, y = 0, width = 0.5, height = 1, just = c("left", "bottom")))
grid.draw(gb_heatmap1)
popViewport()
pushViewport(viewport(x = 0.5, y = 0, width = 0.5, height = 1, just = c("left", "bottom")))
grid.draw(gb_heatmap2)
popViewport()
# refer https://github.com/jokergoo/ComplexHeatmap/issues/110
gb_heatmap1 = grid.grabExpr(draw(h1, padding = unit(c(2, 20, 2, 15), "mm")) )
gb_heatmap2 = grid.grabExpr(draw(h2, padding = unit(c(2, 8, 2, 15), "mm")) )
grid.newpage()
pushViewport(viewport(x = 0, y = 0.95, width = 1, height = 0.48, just = c("left", "top")))
grid.draw(gb_heatmap1)
popViewport()
pushViewport(viewport(x = 0, y = 0, width = 1, height = 0.48, just = c("left", "bottom")))
grid.draw(gb_heatmap2)
popViewport()

lig_tar_heatmap(VISp, 'Glu_Gria2')
lig_tar_heatmap(VISp, 'Glu_Grik5')
lig_tar_heatmap(VISp, 'Cck_Cckbr')

meta_VISp <- data.frame(cell_subclass=VISp@data.signaling$cell_subclass)
rownames(meta_VISp) <- rownames(VISp@data.signaling)
seurat_VISp <- Seurat::CreateSeuratObject(counts = VISp@data, meta.data = meta_VISp)
seurat_VISp$groups <- meta_VISp$cell_subclass
SeuratObject::Idents(seurat_VISp) <- meta_VISp$cell_subclass
levels(seurat_VISp) <- sort(c(sender.names,receiver.names),method='radix')

p1 <- Seurat::VlnPlot(subset(seurat_VISp, idents = sender.names),features = interactionDB$Glu_Gria2$lig_contributor,group.by = 'cell_subclass',pt.size = 0,stack = T,flip = T)
p1$labels$title <- c('Ligand related genes');p1$theme$plot.title$hjust <- 0.5
p2 <- Seurat::VlnPlot(subset(seurat_VISp, idents = receiver.names),features = interactionDB$Glu_Gria2$receptor_subunit,group.by = 'cell_subclass',pt.size = 0) + theme(legend.position = 'none')
p2$labels$title <- paste('Target gene: ',p2$labels$title)
p1 +p2

VlnPlot(subset(seurat_VISp, idents = receiver.names),features = interactionDB$Glu_Gria2$receptor_subunit,group.by = 'cell_subclass',pt.size = 1e-3)

RidgePlot(subset(seurat_VISp, idents = sender.names),features = interactionDB$Glu_Gria2$lig_contributor,group.by = 'cell_subclass')
RidgePlot(subset(seurat_VISp, idents = receiver.names),features = interactionDB$Glu_Gria2$receptor_subunit)

##### Figure2 for ALM  projections
source("~/Documents/CellChat/Rcode_cellchat_VISp_ALM/0220/Projections_ALM.R", encoding = 'UTF-8')
sender.names <- sort(sender.names,method='radix')
receiver.names <- sort(receiver.names,method='radix')
  # set.seed(1);# for weight_count
  # set.seed(1234)
  # df_sender <- readRDS('target_df_ALM.rds')
  # df_receiver <- readRDS('target_df_ALM_proj.rds')
  df_sender <- target_df_ALM
  df_receiver <- target_df_ALM_proj
  df_sender <- df_sender[df_sender$cell_subclass %in% sender.names,]
  df_receiver <- df_receiver[df_receiver$cell_subclass %in% receiver.names,]
  target_df <- rbind(df_sender,df_receiver)
  tmp <- t(as.matrix(target_df[,1:(dim(target_df)[2]-1)]))
  ALM <- createNeuronChat(tmp,DB='mouse', group.by = target_df$cell_subclass);
  ALM <- run_NeuronChat(ALM,sender=sender.names, receiver=receiver.names,N=0,M=100)
  net_aggregated_ALM <- net_aggregation(ALM@net,method='weight_threshold',cut_off = 0.8)
  # net_aggregated_ALM <- net_aggregation(ALM@net[sample(1:length(ALM@net),100)],method='weighted_count',cut_off = 0.8)
  # net_aggregated_ALM <- net_aggregation(ALM@net,method='weighted_count',cut_off = 0.5)
  net.names <- rownames(VISp_projection_mtx)
  net_predicted <- VISp_projection_mtx;net_predicted[sender.names,receiver.names] <- net_aggregated_ALM[sender.names,receiver.names]
  par(mfrow=c(1,2))
  # netVisual_circle_zw(net_aggregated_ALM, title.name = 'Predicted ALM projections', group=group,sources.use =sender.names, targets.use = receiver.names)
  netVisual_circle_neuron(VISp_projection_mtx,group=group,vertex.label.cex=1.5,arrow.size = 0.8,margin = 0.4)
  netVisual_circle_compare(VISp_projection_mtx,(net_predicted[net.names,net.names]/max(net_predicted[net.names,net.names])>0.165)*1,group=group,vertex.label.cex=1.5,arrow.size = 0.8,margin = 0.4)
  library(PRROC)
  PRROC_obj <- roc.curve(scores.class0 = c(net_predicted[sender.names,receiver.names]/max(net_predicted[net.names,net.names])), weights.class0=c(VISp_projection_mtx[sender.names,receiver.names]),curve=TRUE)
  par(mar = c(5, 5, 4, 1),cex.lab=1.5, cex.axis=1.5,cex.main=1.5); plot(PRROC_obj,xlab = 'False positive rate')
  PRROC_obj <- pr.curve(scores.class0 = c(net_predicted[sender.names,receiver.names]/max(net_predicted[net.names,net.names])), weights.class0=c(VISp_projection_mtx[sender.names,receiver.names]),curve=TRUE)
  par(mar = c(5, 5, 4, 1),cex.lab=1.5, cex.axis=1.5,cex.main=1.5);plot(PRROC_obj)

  # PRROC_obj_pr <- pr.curve(scores.class0 = c(net_predicted[sender.names,receiver.names]/max(net_predicted[net.names,net.names])), weights.class0=c(VISp_projection_mtx[sender.names,receiver.names]),curve=TRUE)
  # plot(PRROC_obj_pr)
  accuracy<-c();sensitivity<-c();FPR<-c()
  for(threshold in seq(0,1,1e-3)){
    predicted_net <- (net_predicted/max(net_predicted)>threshold)*1
    accuracy <- c(accuracy,1-sum(abs(predicted_net-VISp_projection_mtx))/length(receiver.names)/length(sender.names))
    sensitivity <- c(sensitivity,sum(c(predicted_net)*(c(VISp_projection_mtx)>0))/sum(c(VISp_projection_mtx)>0))
    FPR <- c(FPR,sum((c(predicted_net[sender.names,receiver.names])*(c(VISp_projection_mtx[sender.names,receiver.names])==0)))/sum((c(VISp_projection_mtx[sender.names,receiver.names])==0)))
  }
  max(accuracy)
  seq(0,1,1e-3)[which(accuracy==max(accuracy))]

plot(FPR,sensitivity,type='l')
VISp_label <- c(VISp_projection_mtx[sender.names,receiver.names])
auroc_all <- rep(0,length(ALM@net))
prroc_all <- rep(0,length(ALM@net))
for(j in 1:length(ALM@net)){
  if(max(ALM@net[[j]])==0){auroc_all[j]=0;prroc_all[j]=0} else{
    VISp_predict <- c(ALM@net[[j]][sender.names,receiver.names])/max(ALM@net[[j]])
    auroc_all[j] <- roc.curve(scores.class0 = VISp_predict, weights.class0=VISp_label,
                              curve=TRUE)$auc
    prroc_all[j] <- pr.curve(scores.class0 = VISp_predict, weights.class0=VISp_label,
                             curve=TRUE)$auc.integral
  }
}
names(auroc_all) <- names(ALM@DB)
auroc_all_sort <- sort(auroc_all,decreasing = T)

weight <- apply(simplify2array(ALM@net),3,sum)
count <- apply(simplify2array(ALM@net),3,FUN=function(x){sum(x>0)})
df_cor <- data.frame(AUROC=auroc_all,AUPRC=prroc_all,Info=weight,count=count);df_cor$interaction <- rownames(df_cor)
df_cor <- df_cor[df_cor$Info>0,]
cor(df_cor$AUROC,df_cor$Info,method = 'spearman')
t <- list()
cor.test(df_cor$AUROC,df_cor$AUPRC,method = 'spearman')
t[[1]] <- cor.test(df_cor$AUROC,df_cor$Info,method = 'spearman')
t[[2]] <- cor.test(df_cor$AUROC,df_cor$count,method = 'spearman')
t[[3]]<- cor.test(df_cor$AUPRC,df_cor$Info,method = 'spearman')
t[[4]] <- cor.test(df_cor$AUPRC,df_cor$count,method = 'spearman')
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
title_rho_p <- sapply(t,FUN=function(x){paste('rho = ',signif(x$estimate,2),'; p = ', signif(x$p.value, 2),sep = '')},simplify = T)
my_lm <- lm(df_cor$AUROC~df_cor$Info)
summary(my_lm)
library(ggplot2)
library(hrbrthemes)
# mtcars dataset is natively available in R
# head(mtcars)
# A basic scatterplot with color depending on Species
p1 <- ggplot(df_cor, aes(x=AUROC, y=Info)) + labs(x ="AUROC", y = "Information flow",title=title_rho_p[1]) +geom_point(size=3)  + theme(text = element_text(size=20),axis.text.x = element_text(angle=0, hjust=1)) + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

p2 <- ggplot(df_cor, aes(x=AUROC, y=count)) + labs(title=title_rho_p[2],
                                                   x ="AUROC", y = "Count of links") +geom_point(size=3)  + theme(text = element_text(size=20),
                                                                                                                  axis.text.x = element_text(angle=0, hjust=1)) +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

p11 <- ggplot(df_cor, aes(x=AUPRC, y=Info)) + labs(title=title_rho_p[3],
                                                   x ="AUPRC", y = "Information flow") +geom_point(size=3)  + theme(text = element_text(size=20),
                                                                                                                    axis.text.x = element_text(angle=0, hjust=1)) +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

p22 <- ggplot(df_cor, aes(x=AUPRC, y=count)) + labs(title=title_rho_p[4],
                                                    x ="AUPRC", y = "Count of links") +geom_point(size=3)  + theme(text = element_text(size=20),
                                                                                                                   axis.text.x = element_text(angle=0, hjust=1)) +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

p1 + p11 +p2 + p22

p3 <- ggplot(df_cor,aes(x = Info,y = reorder(interaction,Info))) + geom_col(fill='orange') + theme(text = element_text(size=18)) + labs(y= "Interaction", x = "Information flow") +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
p4 <- ggplot(df_cor,aes(x = count,y = reorder(interaction,count))) + geom_col(fill='orange') + theme(text = element_text(size=18)) + labs(y= "Interaction", x = "Count of links") +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
p5 <- ggplot(df_cor,aes(x = AUROC,y = reorder(interaction,AUROC))) + geom_col(fill='orange') + theme(text = element_text(size=18)) + labs(y= "Interaction", x = "AUROC") +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) + geom_vline(xintercept=0.5, linetype="dashed", color = 'red', size=0.5)
p6 <- ggplot(df_cor,aes(x = AUPRC,y = reorder(interaction,AUPRC))) + geom_col(fill='orange') + theme(text = element_text(size=18)) + labs(y= "Interaction", x = "AUPRC") +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) + geom_vline(xintercept=sum(VISp_label)/length(sender.names)/length(receiver.names), linetype="dashed", color = 'red', size=0.5)
p5+p6+p3+p4+plot_layout(ncol = 4)



library(ggplot2)
df_auroc <- data.frame(interaction=names(auroc_all_sort),auroc=auroc_all_sort); df_auroc <- df_auroc[df_auroc$auroc!=0,]
ggplot(df_auroc,aes(x = auroc,y = reorder(interaction,auroc))) + geom_col(fill='orange') + theme(text = element_text(size=18)) + labs(y= "Interaction", x = "AUROC")+ geom_vline(xintercept=0.5, linetype="dashed", color = 'red', size=0.5)


## barplot for ligand score & receptor score
col_map=colorspace::diverge_hsv(10)
NeuronChat_object <- ALM
ligand_profile <- NeuronChat_object@ligand.abundance[sender.names,NeuronChat_object@info>0]
ligand_profile <-  ligand_profile[,order(colSums(ligand_profile),decreasing = T)]
target_profile <- NeuronChat_object@target.abundance[receiver.names,NeuronChat_object@info>0]
target_profile <-  target_profile[,order(colSums(target_profile),decreasing = T)]
h1 <- Heatmap(ligand_profile,
              name = "Ligand \n Abundance", cluster_rows = F,cluster_columns = F,column_names_rot = 60,
              row_names_side = 'left',col=col_map,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 12))
h2 <- Heatmap(target_profile,
              name = "Target \n Abundance", cluster_rows = F,cluster_columns = F,column_names_rot = 60,
              row_names_side = 'left',col=col_map,row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 12))
# refer https://github.com/jokergoo/ComplexHeatmap/issues/110
gb_heatmap1 = grid.grabExpr(draw(h1, padding = unit(c(2, 20, 2, 15), "mm")) )
gb_heatmap2 = grid.grabExpr(draw(h2, padding = unit(c(2, 8, 2, 15), "mm")) )
grid.newpage()
pushViewport(viewport(x = 0, y = 0.95, width = 1, height = 0.48, just = c("left", "top")))
grid.draw(gb_heatmap1); popViewport()
pushViewport(viewport(x = 0, y = 0, width = 1, height = 0.48, just = c("left", "bottom")))
grid.draw(gb_heatmap2);popViewport()

lig_tar_heatmap(ALM,'Glu_Gria2')
lig_tar_heatmap(ALM,'Glu_Grik5')
lig_tar_heatmap(ALM,'Nrxn1_Nlgn1')


### supplementary figures related to Figure 2, downsampling VISp
source("~/Documents/CellChat/Rcode_cellchat_VISp_ALM/0220/Projections (1).r", encoding = 'UTF-8')
set.seed(123)
sender.names <- sort(sender.names,method='radix')
receiver.names <- sort(receiver.names,method='radix')
# df_sender <- readRDS('target_df_VISp.rds')
df_sender <- target_df_VISp
df_sender <- df_sender[df_sender$cell_subclass %in% sender.names,]
# df_receiver <- readRDS('target_df_VISp_proj.rds')
df_receiver <- target_df_VISp_proj
df_receiver <- df_receiver[df_receiver$cell_subclass %in% receiver.names,]
target_df <- rbind(df_sender,df_receiver)
tmp <- t(as.matrix(target_df[,1:(dim(target_df)[2]-1)]))
sampling_Rate <- (10:1)/10
VISp_sublist <- list()
roc_list <-  matrix(0,nrow=length(sampling_Rate),ncol=2)
for(j in 1:length(sampling_Rate)){
  sampling_rate <- sampling_Rate[j]
  cell.use.subsampling <- sample(1:dim(tmp)[2],round(dim(tmp)[2]*sampling_Rate[j]))
  VISp_sublist[[j]] <- createNeuronChat(tmp[,cell.use.subsampling],group.by = target_df$cell_subclass[cell.use.subsampling]);
  VISp_sublist[[j]] <- run_NeuronChat(VISp_sublist[[j]],sender=sender.names, receiver=receiver.names,N=0,M=100,fdr=0.05,K=0.5)
  # VISp <- run_NeuronChat(VISp,N=0,M=100)
  net_aggregated_VISp <- net_aggregation(VISp_sublist[[j]]@net,method='weight_threshold',cut_off = 0.8) # weight_threshold
  net.names <- rownames(VISp_projection_mtx)
  net_predicted <- VISp_projection_mtx;net_predicted[sender.names,receiver.names] <- net_aggregated_VISp[sender.names,receiver.names]
  # par(mfrow=c(1,2))
  # # netVisual_circle_zw(net_aggregated_VISp, title.name = 'Predicted VISp projections', group=group,sources.use =sender.names, targets.use = receiver.names,vertex.label.cex=1.5,arrow.size = 0.8)
  # netVisual_circle_zw(VISp_projection_mtx,group=group,vertex.label.cex=1.5,arrow.size = 0.8,margin = 0.4)
  # netVisual_circle_compare(VISp_projection_mtx,(net_predicted[net.names,net.names]/max(net_predicted[net.names,net.names])>0.017)*1,group=group,vertex.label.cex=1.5,arrow.size = 0.8,margin = 0.4)
  library(PRROC)
  AUROC_obj <- roc.curve(scores.class0 = c(net_predicted[sender.names,receiver.names]/max(net_predicted[net.names,net.names])), weights.class0=c(VISp_projection_mtx[sender.names,receiver.names]),curve=TRUE)
  par(mar = c(5, 5, 4, 1),cex.lab=1.5, cex.axis=1.5,cex.main=1.5); plot(AUROC_obj,xlab = 'False positive rate')
  PRROC_obj <- pr.curve(scores.class0 = c(net_predicted[sender.names,receiver.names]/max(net_predicted[net.names,net.names])), weights.class0=c(VISp_projection_mtx[sender.names,receiver.names]),curve=TRUE)
  par(mar = c(5, 5, 4, 1),cex.lab=1.5, cex.axis=1.5,cex.main=1.5);plot(PRROC_obj)
  roc_list[j,1:2]<- c(AUROC_obj$auc,PRROC_obj$auc.integral)
}
df_subsampling <- data.frame(sampling_Rate=sampling_Rate,AUROC=roc_list[,1],AUPRC=roc_list[,2])
plot(sampling_Rate,roc_list[,1],ylim = c(0,1),type = 'p',xlab = 'subsampling rate',ylab = 'AUROC')
g1 <- ggplot(data=df_subsampling, aes(x=sampling_Rate, y=AUROC, group=1)) + xlim(0.05,1) + ylim(0,1) + xlab('subsampling rate')+  theme(text = element_text(size=18)) +
  geom_line(color="red",size=2)+
  geom_point(size=5) #+ geom_hline(yintercept=0.5, linetype="dashed", color = 'gray', size=2)
g2 <- ggplot(data=df_subsampling, aes(x=sampling_Rate, y=AUPRC, group=1)) + xlim(0.05,1) + ylim(0,1) + xlab('subsampling rate')+  theme(text = element_text(size=18)) +
  geom_line(color="red",size=2)+
  geom_point(size=5) #+ geom_hline(yintercept=0.5, linetype="dashed", color = 'gray', size=2)
g1+g2

### supplementary figures related to Figure 2, downsampling VISp
source("~/Documents/CellChat/Rcode_cellchat_VISp_ALM/0220/Projections_ALM.R", encoding = 'UTF-8')
sender.names <- sort(sender.names,method='radix')
receiver.names <- sort(receiver.names,method='radix')
df_sender <- target_df_ALM
df_receiver <- target_df_ALM_proj
# df_sender <- readRDS('target_df_ALM.rds')
df_sender <- df_sender[df_sender$cell_subclass %in% sender.names,]
# df_receiver <- readRDS('target_df_ALM_proj.rds')
df_receiver <- df_receiver[df_receiver$cell_subclass %in% receiver.names,]
target_df <- rbind(df_sender,df_receiver)
tmp <- t(as.matrix(target_df[,1:(dim(target_df)[2]-1)]))
sampling_Rate <- (10:1)/10
ALM_sublist <- list()
roc_list <-  matrix(0,nrow=length(sampling_Rate),ncol=2)
for(j in 1:length(sampling_Rate)){
  sampling_rate <- sampling_Rate[j]
  cell.use.subsampling <- sample(1:dim(tmp)[2],round(dim(tmp)[2]*sampling_Rate[j]))
  ALM_sublist[[j]] <- createNeuronChat(tmp[,cell.use.subsampling],group.by = target_df$cell_subclass[cell.use.subsampling]);
  ALM_sublist[[j]] <- run_NeuronChat(ALM_sublist[[j]],sender=sender.names, receiver=receiver.names,N=0,M=100,fdr=0.05,K=0.5)
  # VISp <- run_NeuronChat(VISp,N=0,M=100)
  net_aggregated_VISp <- net_aggregation(ALM_sublist[[j]]@net,method='weight_threshold',cut_off = 0.8) # weight_threshold
  net.names <- rownames(VISp_projection_mtx)
  net_predicted <- VISp_projection_mtx;net_predicted[sender.names,receiver.names] <- net_aggregated_VISp[sender.names,receiver.names]
  # par(mfrow=c(1,2))
  # # netVisual_circle_zw(net_aggregated_VISp, title.name = 'Predicted VISp projections', group=group,sources.use =sender.names, targets.use = receiver.names,vertex.label.cex=1.5,arrow.size = 0.8)
  # netVisual_circle_zw(VISp_projection_mtx,group=group,vertex.label.cex=1.5,arrow.size = 0.8,margin = 0.4)
  # netVisual_circle_compare(VISp_projection_mtx,(net_predicted[net.names,net.names]/max(net_predicted[net.names,net.names])>0.017)*1,group=group,vertex.label.cex=1.5,arrow.size = 0.8,margin = 0.4)
  library(PRROC)
  AUROC_obj <- roc.curve(scores.class0 = c(net_predicted[sender.names,receiver.names]/max(net_predicted[net.names,net.names])), weights.class0=c(VISp_projection_mtx[sender.names,receiver.names]),curve=TRUE)
  par(mar = c(5, 5, 4, 1),cex.lab=1.5, cex.axis=1.5,cex.main=1.5); plot(AUROC_obj,xlab = 'False positive rate')
  PRROC_obj <- pr.curve(scores.class0 = c(net_predicted[sender.names,receiver.names]/max(net_predicted[net.names,net.names])), weights.class0=c(VISp_projection_mtx[sender.names,receiver.names]),curve=TRUE)
  par(mar = c(5, 5, 4, 1),cex.lab=1.5, cex.axis=1.5,cex.main=1.5);plot(PRROC_obj)
  roc_list[j,1:2]<- c(AUROC_obj$auc,PRROC_obj$auc.integral)
}
df_subsampling_ALM <- data.frame(sampling_Rate=sampling_Rate,AUROC=roc_list[,1],AUPRC=roc_list[,2])
g1 <- ggplot(data=df_subsampling_ALM, aes(x=sampling_Rate, y=AUROC, group=1)) + xlim(0.05,1) + ylim(0,1) + xlab('subsampling rate')+  theme(text = element_text(size=18)) +
  geom_line(color="red",size=2)+
  geom_point(size=5) #+ geom_hline(yintercept=0.5, linetype="dashed", color = 'gray', size=2)
g2 <- ggplot(data=df_subsampling_ALM, aes(x=sampling_Rate, y=AUPRC, group=1)) + xlim(0.05,1) + ylim(0,1) + xlab('subsampling rate')+  theme(text = element_text(size=18)) +
  geom_line(color="red",size=2)+
  geom_point(size=5) #+ geom_hline(yintercept=0.5, linetype="dashed", color = 'gray', size=2)
g1+g2


# ###
# #### Figure 1 test
# cortex_list_tmp <- cortex_list
# LT_list <-c('Glu_Grik4','Glu_Gria1','Cck_Cckbr', 'Glu_Grin2d','Glu_Gria3','Glu_Grm5','Glu_Grm8','Nrxn1_Nlgn1','Glu_Grm3')
# cortex_list_tmp[[1]]@net <- cortex_list_tmp[[1]]@net[LT_list]
# cortex_list_tmp[[2]]@net <- cortex_list_tmp[[2]]@net[LT_list]
# neuronchat_list <- mergeNeuronChat(cortex_list_tmp, add.names = names(cortex_list))
# #### Figure 4c
# g1 <- rankNet_Neuron(neuronchat_list,mode='comparison',measure = c("count"),comparison = 1:2,do.stat = F,tol = 0.1,stacked = F,font.size = 16)
# g2 <- rankNet_Neuron(neuronchat_list,mode='comparison',measure = c("weight"),comparison = 1:2,do.stat = F,tol = 0.1,stacked = F,font.size = 16)#+scale_y_continuous(trans = "log1p",breaks = c(0,1,10,50,200),labels =c(0,1,10,50,200)) #scale_y_continuous(breaks = log2(c(1,10,120)), labels =log2(c(1,10,120))) ## + scale_y_log10()
# g1+g2
