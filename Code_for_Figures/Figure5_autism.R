library(NeuronChat)
require(Seurat)
require(data.table)
library(Matrix)
library(nichenetr)
library(circlize)
set.seed(1234)
setwd('/Users/weizhao/Documents/CellChat/ASD/')
mat <- Matrix::readMM(file='/Users/weizhao/Documents/CellChat/ASD/rawMatrix/matrix.mtx')
barcode.names <- read.delim('/Users/weizhao/Documents/CellChat/ASD/rawMatrix/barcodes.tsv',header = FALSE,stringsAsFactors = FALSE)
feature.names <- read.delim('/Users/weizhao/Documents/CellChat/ASD/rawMatrix/genes.tsv',header = FALSE,stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2
meta <- read.table("/Users/weizhao/Documents/CellChat/ASD/rawMatrix/meta.txt", header=T, sep="\t", as.is=T, row.names=1)

library.size <- Matrix::colSums(mat)
mat_tmp <- t(t(mat)/ library.size) * 10000
mat_norm <- log1p(mat_tmp)

setwd('/Users/weizhao/Documents/CellChat/Rcode_cellchat_VISp_ALM/')
interactionDB <- readRDS('interactionDB_human.rds')
interactionDB <- interactionDB[!duplicated(names(interactionDB))]
interactionDB <- interactionDB[sort(names(interactionDB),method='radix')]
target_genes <- c()
for(j in 1:length(interactionDB)){
  target_genes <- c(target_genes,interactionDB[[j]]$lig_contributor,interactionDB[[j]]$receptor_subunit,interactionDB[[j]]$targets_up,interactionDB[[j]]$targets_down,interactionDB[[j]]$activator,interactionDB[[j]]$inhibitor, interactionDB[[j]]$interactors)
}
target_genes <- unique(target_genes)
target_mat <- mat_norm[which(rownames(mat_norm) %in% target_genes),]

target_df_ASD_Control <- as.data.frame(as.matrix(t(target_mat)));colnames(target_df_ASD_Control) <-rownames(mat_norm)[which(rownames(mat_norm) %in% target_genes)]
target_df_ASD_Control$cell_subclass <- meta$cluster
saveRDS(target_df_ASD_Control,file='autism_mat.rds')
saveRDS(meta,file='autism_meta.rds')

##
target_df_ASD_Control <- readRDS('autism_mat.rds');
meta <- readRDS('autism_meta.rds');
##### Figure5
#############
set.seed(1234)
ASD_Control <- c('ASD','Control')
ASD_list <- lapply(ASD_Control, function(y){
  region_name <- y;
  cell_idx <- which(meta$diagnosis %in% region_name); cell_idx <-sample(cell_idx,1e4)
  target_df_single  <- target_df_ASD_Control[cell_idx,]
  meta_tmp <- meta[cell_idx,];
  x <- createNeuronChat( t(as.matrix(target_df_single[,1:(dim(target_df_single)[2]-1)])),DB='human',group.by = target_df_single$cell_subclass,meta=meta_tmp);
  x <- run_NeuronChat(x,M=100)
  net_aggregated_x <- net_aggregation(x@net,method='weight')
  netVisual_circle_neuron(net_aggregated_x, title.name = y)
  return(x)
})
names(ASD_list) <- ASD_Control
##### Figure5b
#############
net_aggregated_x <- net_aggregation(ASD_list[[1]]@net,'weight')
net_aggregated_y <- net_aggregation(ASD_list[[2]]@net,'weight')
par(mfrow=c(1,2),mar=c(2,2,2,4))
netVisual_circle_neuron((abs(net_aggregated_x-net_aggregated_y)+(net_aggregated_x-net_aggregated_y))/2,title.name = 'Upregulated communications in ASD',margin = 0.25,vertex.label.cex = 1.5)
netVisual_circle_neuron((abs(net_aggregated_x-net_aggregated_y)-(net_aggregated_x-net_aggregated_y))/2,title.name = 'Downregulated communications in ASD',margin = 0.25,vertex.label.cex = 1.5)

##### Figure5e
#############
neuronchat_list <- mergeNeuronChat(ASD_list, add.names = names(ASD_list))
neuronchat_list <- computeNetSimilarityPairwise_Neuron(neuronchat_list, slot.name = "net", type = "functional",comparison = NULL)
#> Compute signaling network similarity for datasets 1 2
neuronchat_list <- netEmbedding(neuronchat_list,slot.name = "net_analysis", type = "functional",comparison =NULL)
#> Manifold learning of the signaling networks for datasets 1 2
neuronchat_list <- netClustering(neuronchat_list, slot.name = "net_analysis", type = "functional",comparison = NULL)
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise_Neuron(neuronchat_list, slot.name = "net_analysis", type = "functional", label.size = 3.5,comparison=NULL,pathway.remove.show = FALSE,pathway.labeled=c())
netVisual_embeddingPairwiseZoomIn_Neuron(neuronchat_list, slot.name = "net_analysis", type = "functional", label.size = 3.5,comparison=NULL,nCol=2)

##### Figure5d
#############
gg1 <- rankNet_Neuron(neuronchat_list,mode='comparison',measure = c("count"),comparison = c(1,2),do.stat = F,font.size = 12)
gg2 <- rankNet_Neuron(neuronchat_list,mode='comparison',measure = c("weight"),comparison = c(1,2),do.stat = F,font.size = 12)
gg1 + gg2

##### Figure5a
#############
gg1 <- compareInteractions_Neuron(neuronchat_list,measure = c("count"),comparison = c(1,2),show.legend = F, group=c(1,2),size.text = 20)
gg2 <- compareInteractions_Neuron(neuronchat_list, measure = c("weight"), comparison = c(1,2), show.legend = F, group = c(1,2),size.text = 20)
gg1 + gg2


# cluster_info <- neuronchat_list@net_analysis$similarity$functional$group$`1-2`
# cluster_info_region <- cluster_info;names(cluster_info_region) <- sub(".*--", "", names(cluster_info_region))
# cluster_info_pathway <- cluster_info;names(cluster_info_pathway) <- sub("--.*","", names(cluster_info_pathway))
# df_info <- data.frame(cluster=cluster_info,region=names(cluster_info_region),pathway=names(cluster_info_pathway))
# df_info <- df_info[order(df_info$pathway),]

##### Figure5c
#############
info <- data.frame(control=ASD_list[[2]]@info,ASD=ASD_list[[2]]@info)
interactionDB <- ASD_list[[1]]@DB
info$ratio <- info$ASD/info$control;info$interaction_names <- names(interactionDB)
info <- info[!((info$control^2+info$ASD^2)==0),]
info$control_percent <- info$control/(info$control+info$ASD)
info$ASD_percent <- info$ASD/(info$control+info$ASD)
info <- info[order(info$ASD_percent,decreasing = T),]
info<- info[!duplicated(info$interaction_names),]

# info_long <- data.frame(value=c(info$control,info$ASD),interaction=c(info$interaction_names,info$interaction_names),percent = c(info$control_percent,info$ASD_percent),condition=kronecker(1:2, rep(1,dim(info)[1])))
# info_long$diagnosis <- c('Control','ASD')[info_long$condition]
# info_long$interaction <- factor(info_long$interaction,levels=info$interaction)
# ggplot(info_long, aes(fill=diagnosis, y=interaction,x=percent)) +
#   geom_bar(position="fill", stat="identity") + labs(y= "Interaction", x = "Relative information flow") +  theme(text = element_text(size=20),
#                                                                                                                 axis.text.y = element_text(angle=0, hjust=1))
# ggplot(info_long, aes(fill=diagnosis, y=interaction,x=value)) +
#   geom_bar(position="dodge", stat="identity") + labs(y= "Interaction", x = "Relative information flow") +  theme(text = element_text(size=20),
#                                                                                                                 axis.text.y = element_text(angle=0, hjust=1))
interaction_use <- info$interaction
interaction_use_idx <- match(interaction_use,names(interactionDB))
net.rownames <- rownames(ASD_list[[1]]@net[[1]])
x <- ASD_list[[2]]
y <- ASD_list[[1]]
cell_sending_x <- matrix(0,nrow=length(net.rownames),ncol=length(interaction_use_idx));rownames(cell_sending_x) <- net.rownames;colnames(cell_sending_x) <- interaction_use
cell_receiving_x <- cell_sending_x
for(j in 1:length(interaction_use_idx)){
  cell_sending_x[1:length(net.rownames),j] <- rowSums(x@net[[interaction_use_idx[j]]])
  cell_receiving_x[1:length(net.rownames),j] <- colSums(x@net[[interaction_use_idx[j]]])
}
cell_sending_y<- cell_sending_x
cell_receiving_y <- cell_sending_y
for(j in 1:length(interaction_use_idx)){
  cell_sending_y[1:length(net.rownames),j] <- rowSums(y@net[[interaction_use_idx[j]]])
  cell_receiving_y[1:length(net.rownames),j] <- colSums(y@net[[interaction_use_idx[j]]])
}
diff_sending <- (cell_sending_y-cell_sending_x)#/(cell_sending_y+cell_sending_x+1e-16)
colmap <- circlize::colorRamp2(c(-0.1,0,0.5), c("blue","slategray1","red"), transparency = 0, space = "LAB")
h1 <- Heatmap(diff_sending,name = "Diff",cluster_rows = FALSE,cluster_columns=FALSE,row_names_side='left',col=colmap,column_names_rot = 50)
diff_receiving <- (cell_receiving_y-cell_receiving_x)#/(cell_receiving_y+cell_receiving_x+1e-16)
colmap <- circlize::colorRamp2(c(-0.1,0,0.5), c("blue","slategray1","red"), transparency = 0, space = "LAB")
h2 <- Heatmap(diff_receiving,name = "Diff",cluster_rows = FALSE,cluster_columns=FALSE,row_names_side='left',col=colmap,column_names_rot = 50)
draw(h1, column_title = 'Differential outgoing signal strength', column_title_gp = gpar(fontsize = 16))
draw(h2, column_title = 'Differential incoming signal strength', column_title_gp = gpar(fontsize = 16))
gb_h1 <- grid.grabExpr(draw(h1,column_title = 'Differential outgoing signal strength', column_title_gp = gpar(fontsize = 16)))
gb_h2 <- grid.grabExpr(draw(h2,column_title = 'Differential incoming signal strength', column_title_gp = gpar(fontsize = 16)))
grid.newpage()
pushViewport(viewport(x = 0, y = 1, width = 0.95, height =0.5, just = c("left", "top")))
grid.draw(gb_h1);popViewport()
pushViewport(viewport(x = 0, y = 0.5, width = 0.95, height = 0.5,just = c("left", "top")))
grid.draw(gb_h2);popViewport()



