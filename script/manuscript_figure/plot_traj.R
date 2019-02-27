library(slingshot)
library(ggplot2)
library(RColorBrewer)
library(CellBench)
library(scran)
library(scater)
library(mclust)
library(monocle)
library(dplyr)

# feature selection
scran_high_var = function(sce,topn=1000){
  var.fit <- trendVar(sce, method="loess", use.spikes=FALSE)
  var.out <- decomposeVar(sce, var.fit)
  hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
  rowData(sce)$hi_var = FALSE
  rowData(sce)$hi_var[rownames(rowData(sce)) %in% rownames(hvg.out)] = TRUE
  return(sce)
}

prep_RNA_traj_order = function(sce){
  #sce$group = paste(colData(sce)$H2228_prop,colData(sce)$H1975_prop,colData(sce)$HCC827_prop,sep="_")
  sce$H2228_to_H1975 = NA
  sce$H2228_to_H1975[sce$group=="1 0 0"] = 0
  sce$H2228_to_H1975[sce$group=="0.68 0.16 0.16"] = 1
  sce$H2228_to_H1975[sce$group=="0.33 0.33 0.33"] = 2
  sce$H2228_to_H1975[sce$group=="0.16 0.68 0.16"] = 3
  sce$H2228_to_H1975[sce$group=="0 1 0"] = 4
  
  sce$H2228_to_HCC827 = NA
  sce$H2228_to_HCC827[sce$group=="1 0 0"] = 0
  sce$H2228_to_HCC827[sce$group=="0.68 0.16 0.16"] = 1
  sce$H2228_to_HCC827[sce$group=="0.33 0.33 0.33"] = 2
  sce$H2228_to_HCC827[sce$group=="0.16 0.16 0.68"] = 3
  sce$H2228_to_HCC827[sce$group=="0 0 1"] = 4
  
  sce$H1975_to_HCC827 = NA
  sce$H1975_to_HCC827[sce$group=="0 1 0"] = 0
  sce$H1975_to_HCC827[sce$group=="0.16 0.68 0.16"] = 1
  sce$H1975_to_HCC827[sce$group=="0.33 0.33 0.33"] = 2
  sce$H1975_to_HCC827[sce$group=="0.16 0.16 0.68"] = 3
  sce$H1975_to_HCC827[sce$group=="0 0 1"] = 4
  return(sce)
}

prep_traj_order = function(sce){
  #sce$group = paste(colData(sce)$H2228,colData(sce)$H1975,colData(sce)$HCC827,sep="_")
  sce$H2228_to_H1975 = NA
  sce$H2228_to_H1975[sce$group=="9 0 0"] = 0
  sce$H2228_to_H1975[sce$group=="7 1 1"] = 1
  sce$H2228_to_H1975[sce$group=="5 2 2"] = 2
  sce$H2228_to_H1975[sce$group=="3 3 3"] = 3
  sce$H2228_to_H1975[sce$group=="2 5 2"] = 4
  sce$H2228_to_H1975[sce$group=="1 7 1"] = 5
  sce$H2228_to_H1975[sce$group=="0 9 0"] = 6
  
  sce$H2228_to_HCC827 = NA
  sce$H2228_to_HCC827[sce$group=="9 0 0"] = 0
  sce$H2228_to_HCC827[sce$group=="7 1 1"] = 1
  sce$H2228_to_HCC827[sce$group=="5 2 2"] = 2
  sce$H2228_to_HCC827[sce$group=="3 3 3"] = 3
  sce$H2228_to_HCC827[sce$group=="2 2 5"] = 4
  sce$H2228_to_HCC827[sce$group=="1 1 7"] = 5
  sce$H2228_to_HCC827[sce$group=="0 0 9"] = 6
  
  sce$H1975_to_HCC827 = NA
  sce$H1975_to_HCC827[sce$group=="0 9 0"] = 0
  sce$H1975_to_HCC827[sce$group=="1 7 1"] = 1
  sce$H1975_to_HCC827[sce$group=="2 5 2"] = 2
  sce$H1975_to_HCC827[sce$group=="3 3 3"] = 3
  sce$H1975_to_HCC827[sce$group=="2 2 5"] = 4
  sce$H1975_to_HCC827[sce$group=="1 1 7"] = 5
  sce$H1975_to_HCC827[sce$group=="0 0 9"] = 6
  
  sce = sce[,colData(sce)$traj =="YES"] # there are cells mixed with two cell lines and are not in the trajectory path
  return(sce)
}


#RNAmix_after_trajectory <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/RNAmix_after_trajectory.Rds")
RNAmix_all_after_imputation <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/RNAmix_all_after_imputation.Rds")
cellmix_all_after_imputation <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/cellmix_all_after_imputation.Rds")

#### slingshot
start_grp="1 0 0"
sce = RNAmix_all_after_imputation[RNAmix_all_after_imputation$data=="RNAmix_CELseq2"& RNAmix_all_after_imputation$norm_method== "Linnorm" & RNAmix_all_after_imputation$impute_method== "SAVER",]$result[[1]]
print(assayNames(sce))
sce = prep_RNA_traj_order(sce)
sce_de_traj = runPCA(sce)
#kmeans_clu = kmeans(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"),centers=num_k,iter.max = 10000)
cl1 <- Mclust(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"))$classification
colData(sce_de_traj)$GMM_cluster= as.factor(cl1)
tmp = table(colData(sce_de_traj)$GMM_cluster[colData(sce_de_traj)$group==start_grp])  # specify H2228 as root state.
tmp = tmp[order(tmp,decreasing = T)]

slingshot_lin <- slingshot::getLineages(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"), colData(sce_de_traj)$GMM_cluster,start.clus = names(tmp)[1])
slingshot_crv <- slingshot::getCurves(slingshot_lin)
slingshot_pseudo <- slingshot::slingPseudotime(slingshot_crv)
print(slingshot_pseudo)
pdf("slingshot_Linnorm_SAVER_RNAmix.pdf")
print(plot(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"), 
           col = brewer.pal(9,"YlGnBu")[3:8][as.factor(sce_de_traj$H2228_prop)], 
           asp = 1, pch = 16,xaxt = "n",yaxt = "n"))
print(lines(slingshot_crv, lwd = 3))
print(pairs(slingshot_pseudo))
#ggplot(data=NULL,aes(x=reducedDim(sce_de_traj,"PCA")[,1],y=reducedDim(sce_de_traj,"PCA")[,2],col=slingshot_crv[,1]))+geom_point()
#ggplot(data=NULL,aes(x=reducedDim(sce_de_traj,"PCA")[,1],y=reducedDim(sce_de_traj,"PCA")[,2],col=slingshot_crv[,2]))+geom_point()
dev.off()


start_grp="9 0 0"
sce = cellmix_all_after_imputation[cellmix_all_after_imputation$data=="cellmix3"& cellmix_all_after_imputation$norm_method== "scone" & cellmix_all_after_imputation$impute_method== "SAVER",]$result[[1]]
print(assayNames(sce))
sce = prep_traj_order(sce)
sce_de_traj = runPCA(sce)
#kmeans_clu = kmeans(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"),centers=num_k,iter.max = 10000)
cl1 <- Mclust(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"))$classification
colData(sce_de_traj)$GMM_cluster= as.factor(cl1)
tmp = table(colData(sce_de_traj)$GMM_cluster[colData(sce_de_traj)$group==start_grp])  # specify H2228 as root state.
tmp = tmp[order(tmp,decreasing = T)]

slingshot_lin <- slingshot::getLineages(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"), colData(sce_de_traj)$GMM_cluster,start.clus = names(tmp)[1])
slingshot_crv <- slingshot::getCurves(slingshot_lin)
slingshot_pseudo <- slingshot::slingPseudotime(slingshot_crv)
print(slingshot_pseudo)
pdf("slingshot_scone_SAVER_cellmix3.pdf")
print(plot(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"), 
           col = brewer.pal(9,"YlGnBu")[3:8][as.factor(sce_de_traj$H2228)], 
           asp = 1, pch = 16,xaxt = "n",yaxt = "n"))
print(lines(slingshot_crv, lwd = 3))
#ggplot(data=NULL,aes(x=reducedDim(sce_de_traj,"PCA")[,1],y=reducedDim(sce_de_traj,"PCA")[,2],col=slingshot_crv[,1]))+geom_point()
#ggplot(data=NULL,aes(x=reducedDim(sce_de_traj,"PCA")[,1],y=reducedDim(sce_de_traj,"PCA")[,2],col=slingshot_crv[,2]))+geom_point()
dev.off()

###### Monocle2
sce = RNAmix_all_after_imputation[RNAmix_all_after_imputation$data=="RNAmix_CELseq2"& RNAmix_all_after_imputation$norm_method== "scran" & RNAmix_all_after_imputation$impute_method== "DrImpute",]$result[[1]]
start_grp="1 0 0"
sce = prep_RNA_traj_order(sce)
sce = scran_high_var(sce)
pd = new("AnnotatedDataFrame", data = as.data.frame(colData(sce)))
cds = new("CellDataSet", exprs = logcounts(sce), phenoData = pd, expressionFamily = VGAM::negbinomial.size())
sizeFactors(cds)=1 # assume the data has been normalized.
if("hi_var" %in% colnames(rowData(sce))){
  cds = setOrderingFilter(cds, ordering_genes = rownames(sce)[rowData(sce)$hi_var])
}
cds = reduceDimension(cds, method = "DDRTree",norm_method="none",pseudo_expr=0)
cds = orderCells(cds)
tmp = table(pData(cds)$State[pData(cds)$group==start_grp])  # specify H2228 as root state.
tmp = tmp[order(tmp,decreasing = T)]
cds = orderCells(cds,root_state=names(tmp)[1]) # reorder the cells


reduced_dim_coords <- as.data.frame(t(monocle::reducedDimS(cds)))
colnames(reduced_dim_coords) = c("Dim1","Dim2")

ica_space_df <- Matrix::t(monocle::reducedDimK(cds)) %>%
  as.data.frame() %>%
  mutate(sample_name = rownames(.), sample_state = rownames(.))

colnames(ica_space_df) = c("prin_graph_dim_1","prin_graph_dim_2","sample_name","sample_state")
dp_mst <- minSpanningTree(cds)
edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  select_(source = "from", target = "to") %>%
  left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
  left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")


p_monocle1 = ggplot()+
  geom_point(data=reduced_dim_coords,aes(x=Dim1,y=Dim2,col=factor(sce$H2228_prop)),show.legend = F,size=3,alpha=0.9)+
  scale_color_manual(values=brewer.pal(9,"YlGnBu")[3:8])+
  geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), linetype="solid", na.rm=TRUE, data=edge_df)+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
pdf("Monocle2_scran_DrImpute_RNAmix.pdf")
p_monocle1
dev.off()
p_monocle1




start_grp="9 0 0"
sce = cellmix_all_after_imputation[cellmix_all_after_imputation$data=="cellmix3"& cellmix_all_after_imputation$norm_method== "scran" & cellmix_all_after_imputation$impute_method== "no_impute",]$result[[1]]

sce = prep_traj_order(sce)
sce = scran_high_var(sce)
pd = new("AnnotatedDataFrame", data = as.data.frame(colData(sce)))
cds = new("CellDataSet", exprs = logcounts(sce), phenoData = pd, expressionFamily = VGAM::negbinomial.size())
sizeFactors(cds)=1 # assume the data has been normalized.
if("hi_var" %in% colnames(rowData(sce))){
  cds = setOrderingFilter(cds, ordering_genes = rownames(sce)[rowData(sce)$hi_var])
}
cds = reduceDimension(cds, method = "DDRTree",norm_method="none",pseudo_expr=0)
cds = orderCells(cds)
tmp = table(pData(cds)$State[pData(cds)$group==start_grp])  # specify H2228 as root state.
tmp = tmp[order(tmp,decreasing = T)]
cds = orderCells(cds,root_state=names(tmp)[1]) # reorder the cells


reduced_dim_coords <- as.data.frame(t(monocle::reducedDimS(cds)))
colnames(reduced_dim_coords) = c("Dim1","Dim2")

ica_space_df <- Matrix::t(monocle::reducedDimK(cds)) %>%
  as.data.frame() %>%
  mutate(sample_name = rownames(.), sample_state = rownames(.))

colnames(ica_space_df) = c("prin_graph_dim_1","prin_graph_dim_2","sample_name","sample_state")
dp_mst <- minSpanningTree(cds)
edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  select_(source = "from", target = "to") %>%
  left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
  left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")


p_monocle2 = ggplot()+
  geom_point(data=reduced_dim_coords,aes(x=Dim1,y=Dim2,col=factor(sce$H2228)),show.legend = F,size=3,alpha=0.9)+
  scale_color_manual(values=brewer.pal(9,"YlGnBu")[2:8])+
  geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), linetype="solid", na.rm=TRUE, data=edge_df)+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
pdf("Monocle2_scran_cellmix3.pdf")
p_monocle2
dev.off()
p_monocle2



#### SLICER

library(SLICER)
library(lle)
sce = RNAmix_all_after_imputation[RNAmix_all_after_imputation$data=="RNAmix_CELseq2"& RNAmix_all_after_imputation$norm_method== "Linnorm" & RNAmix_all_after_imputation$impute_method== "DrImpute",]$result[[1]]
start_grp="1 0 0"
sce = prep_RNA_traj_order(sce)
sce = scran_high_var(sce)

#genes = select_genes(traj)
k = select_k(t(logcounts(sce[rowData(sce)$hi_var,])), kmin=2)
traj_lle = lle(t(logcounts(sce[rowData(sce)$hi_var,])), m=2, k=k)$Y
traj_graph = conn_knn_graph(traj_lle,5)
#ends = find_extreme_cells(traj_graph, traj_lle)
#start = 100
#cells_ordered = cell_order(traj_graph, start)
#branches = assign_branches(traj_graph,start)

pdf("SLICER_Linnorm_DrImpute_RNAmix.pdf")
ggplot(data=NULL,aes(x=traj_lle[,1],y=traj_lle[,2],col=factor(sce$H2228_prop)))+
  geom_point(show.legend = F,size=3,alpha=0.9)+
  scale_color_manual(values=brewer.pal(9,"YlGnBu")[3:8])+
  labs(x="Dim1",y="Dim2")+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()



start_grp="9 0 0"
sce = cellmix_all_after_imputation[cellmix_all_after_imputation$data=="cellmix3"& cellmix_all_after_imputation$norm_method== "BASiCS" & cellmix_all_after_imputation$impute_method== "SAVER",]$result[[1]]

sce = prep_traj_order(sce)
sce = scran_high_var(sce)

k = select_k(t(logcounts(sce[rowData(sce)$hi_var,])), kmin=2)
traj_lle = lle(t(logcounts(sce[rowData(sce)$hi_var,])), m=2, k=k)$Y
#traj_graph = conn_knn_graph(traj_lle,5)
#ends = find_extreme_cells(traj_graph, traj_lle)
#start_grp="9_0_0"
#H2228_st = which(colData(sce_SC2_qc_traj)$group == start_grp)[1]
#cells_ordered = cell_order(traj_graph, H2228_st)
#branches = assign_branches(traj_graph,H2228_st)

pdf("SLICER_BASiCS_SAVER_cellmix3.pdf")
ggplot(data=NULL,aes(x=traj_lle[,1],y=traj_lle[,2],col=factor(sce$H2228)))+
  geom_point(show.legend = F,size=3)+
  scale_color_manual(values=brewer.pal(9,"YlGnBu")[2:8])+
  labs(x="Dim1",y="Dim2")+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()