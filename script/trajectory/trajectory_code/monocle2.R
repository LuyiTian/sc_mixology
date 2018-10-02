# monocle2
library(monocle)
method="monocle-DDRTree"
get_max_score = function(col_anno, states){
  col_st = col_anno[col_anno$State %in% states,]
  
  cor1 = cor(col_st$Pseudotime, col_st$H2228_to_H1975, use = "pairwise.complete.obs")
  cor2 = cor(col_st$Pseudotime, col_st$H2228_to_HCC827, use = "pairwise.complete.obs")
  cor3 = cor(col_st$Pseudotime, col_st$H1975_to_HCC827, use = "pairwise.complete.obs")
  
  ov1 = sum(col_anno[!is.na(col_anno$H2228_to_H1975),"State"] %in% states)/sum(!is.na(col_anno$H2228_to_H1975))
  ov2 = sum(col_anno[!is.na(col_anno$H2228_to_HCC827),"State"] %in% states)/sum(!is.na(col_anno$H2228_to_HCC827))
  ov3 = sum(col_anno[!is.na(col_anno$H1975_to_HCC827),"State"] %in% states)/sum(!is.na(col_anno$H1975_to_HCC827))
  
  res_df = data.frame(corr=abs(c(cor1,cor2,cor3)),overlap=c(ov1,ov2,ov3))
  #print(res_df)
  res_df = res_df[order(res_df$overlap,decreasing = T),]
  return(res_df[1,])
}


monocle2_sup_order = function(sce, dr_method="DDRTree", sup_genes=c(),start_grp="9_0_0"){
  pd <- new("AnnotatedDataFrame", data = as.data.frame(colData(sce)))
  cds <- newCellDataSet(counts(sce), phenoData = pd,expressionFamily=negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  if(length(sup_genes)>0){
    print(paste("supervised with known gene list:",length(sup_genes)))
    cds = setOrderingFilter(cds, ordering_genes = sup_genes)
    cds = reduceDimension(cds, method = dr_method)
    cds = orderCells(cds)
  }else{
    print("unsupervised. clustering and get DE genes.")
    disp_table = dispersionTable(cds)
    unsup_clustering_genes = subset(disp_table, mean_expression >= 0.1)
    cds = setOrderingFilter(cds, unsup_clustering_genes$gene_id)
    cds = reduceDimension(cds, max_components = 2, num_dim = 4, # the variations on benchmark data are mostly in first few dims.
                            reduction_method = 'tSNE', verbose = F)
    cds = clusterCells(cds, num_clusters = 10)
    clustering_DEG_genes = differentialGeneTest(cds,
                           fullModelFormulaStr = '~Cluster',
                           cores = 4)
    unsup_clustering_genes = row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
    cds = setOrderingFilter(cds,
                        ordering_genes = unsup_clustering_genes)
    print(plot_cell_clusters(cds, color_by = 'as.factor(Cluster)'))
    cds = reduceDimension(cds, method = dr_method)
    cds = orderCells(cds)
  }
  tmp = table(pData(cds)$State[pData(cds)$group==start_grp])  # specify H2228 as root state.
  tmp = tmp[order(tmp,decreasing = T)]
  cds = orderCells(cds,root_state=names(tmp)[1]) # reorder the cells
  max_score_df = list()
  for(i in 1:length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points)){
    cds_reduced = buildBranchCellDataSet(cds, branch_point=i) # get branch info
    col_data_df = as.data.frame(pData(cds))
    max_score_df[[i]] = lapply(1:length(unique(pData(cds_reduced)$Branch)),function(j){
      state1 = table(pData(cds_reduced)$State[pData(cds_reduced)$Branch == unique(pData(cds_reduced)$Branch)[j]])
      state1 = state1[state1>1]
      print(state1)
      max_score_df1 = get_max_score(col_data_df, names(state1))
    })
    max_score_df[[i]] = Reduce(rbind,max_score_df[[i]])
    print(max_score_df[[i]])
  }
  tmp = unlist(lapply(max_score_df,function(x){colMeans(x)[2]}))
  #print(max_score_df[[which(tmp==max(tmp))]])
  return(max_score_df[[which(tmp==max(tmp))]])
}

ptm <- proc.time()

load("~/Dropbox/research/benchmark/rdata/9cellmix_qc.RData")
source('~/Dropbox/research/benchmark/trajectory_code/util_func.R')

get_traj_rest = function(sce,start_grp="9_0_0"){
  if(start_grp=="1_0_0"){
    sce = prep_RNA_traj_order(sce)
    sce = scran_norm(sce)
    sce_traj = sce
  }
  else{
    sce = prep_traj_order(sce)
    sce = scran_norm(sce)
    sce_traj = sce[,colData(sce)$traj =="YES"]
  }

  max_score = list()
  high_var_genes = scran_high_var(sce_traj)
  for(i in 1:10){
    set.seed(1000*i) # so the result is reproducible
    pr = (1000:1)+100
    pr = pr/sum(pr)
    high_var_sel = sample(high_var_genes,size=500,prob=pr)
    max_score[[i]] = monocle2_sup_order(sce_traj, sup_genes=high_var_sel)
  }
  max_score = Reduce(rbind,max_score)
  return(max_score)
}

max_score_SC1 = get_traj_rest(sce_SC1_qc)
max_score_SC1$method=method
max_score_SC1$design="cellmix"
max_score_SC1$dataset="cellmix1"

max_score_SC2 = get_traj_rest(sce_SC2_qc)
max_score_SC2$method=method
max_score_SC2$design="cellmix"
max_score_SC2$dataset="cellmix2"

max_score_SC3 = get_traj_rest(sce_SC3_qc)
max_score_SC3$method=method
max_score_SC3$design="cellmix"
max_score_SC3$dataset="cellmix3"

max_score_SC4 = get_traj_rest(sce_SC4_qc)
max_score_SC4$method=method
max_score_SC4$design="cellmix"
max_score_SC4$dataset="cellmix4"


load("~/Dropbox/research/benchmark/rdata/mRNAmix_qc.RData")

max_score_sce2 = get_traj_rest(sce2_qc,start_grp="1_0_0")
max_score_sce2$method=method
max_score_sce2$design="RNAmix"
max_score_sce2$dataset="RNAmix_CEL-seq2"

max_score_sce8 = get_traj_rest(sce8_qc,start_grp="1_0_0")
max_score_sce8$method=method
max_score_sce8$design="RNAmix"
max_score_sce8$dataset="RNAmix_Sort-seq"


final_res = rbind(max_score_SC1,max_score_SC2,max_score_SC3,max_score_SC4,
                  max_score_sce2, max_score_sce8)

write.csv(final_res,file=paste0("~/Dropbox/research/benchmark/traj_result/trajectory_result_",method,".csv"),row.names = F)


final_time = (proc.time() - ptm)[3]

write.csv(data.frame(method=method,time_taken=final_time),file=paste0("~/Dropbox/research/benchmark/traj_result/trajectory_time_",method,"_.csv"),row.names = F)

  