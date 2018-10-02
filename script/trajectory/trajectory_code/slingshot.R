# slingshot
method="slingshot"
library(slingshot)
library(scran)
library(scater)

get_max_score = function(col_anno, pseudo){
  col_st = col_anno
  
  cor1 = cor(pseudo, col_st$H2228_to_H1975, use = "pairwise.complete.obs")
  cor2 = cor(pseudo, col_st$H2228_to_HCC827, use = "pairwise.complete.obs")
  
  ov1 = sum(!is.na(pseudo[!is.na(col_anno$H2228_to_H1975)]))/sum(!is.na(col_anno$H2228_to_H1975))
  ov2 = sum(!is.na(pseudo[!is.na(col_anno$H2228_to_HCC827)]))/sum(!is.na(col_anno$H2228_to_HCC827))
  
  res_df = data.frame(corr=abs(c(cor1,cor2)),overlap=c(ov1,ov2))
  print(res_df)
  res_df = res_df[order(res_df$overlap,decreasing = T),]
  return(res_df[1,])
}



slingshot_sup_order = function(sce, start_grp="9_0_0",num_k=10){
  res_all=list()
  high_var_genes = scran_high_var(sce)
  for(i in 1:10){
    set.seed(1000*i) # so the result is reproducible
    pr = (1000:1)+100
    pr = pr/sum(pr)
    high_var_sel = sample(high_var_genes,size=500,prob=pr)
    sce_de_traj = plotPCA(sce,return_SCE=T, draw_plot=T, rerun=T,run_args=list(feature_set=high_var_sel))
    kmeans_clu = kmeans(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"),centers=num_k,iter.max = 10000)
    colData(sce_de_traj)$kmeans_cluster= as.factor(kmeans_clu$cluster)
    tmp = table(colData(sce_de_traj)$kmeans_cluster[colData(sce_de_traj)$group==start_grp])  # specify H2228 as root state.
    tmp = tmp[order(tmp,decreasing = T)]
  
    slingshot_lin <- getLineages(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"), colData(sce_de_traj)$kmeans_cluster,start.clus = names(tmp)[1])
    slingshot_crv <- getCurves(slingshot_lin)
    slingshot_pseudo <- pseudotime(slingshot_crv)
    col_data_df = as.data.frame(colData(sce_de_traj))
    score_df = apply(slingshot_pseudo,2,function(x){get_max_score(col_data_df, x)})
    #print(plot(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"), col = brewer.pal(num_k,"Set3")[colData(sce_de_traj)$kmeans_cluster], asp = 1, pch = 16))
    #print(lines(slingshot_lin, lwd = 3,show.constraints = TRUE))
    #print(plot(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"), 
    #           col = brewer.pal(num_k,"Set3")[colData(sce_de_traj)$kmeans_cluster], 
    #           asp = 1, pch = 16))
    score_df = Reduce(rbind,score_df)
    score_df = score_df[order(score_df$corr,decreasing = T),]
    #print(lines(slingshot_crv, lwd = 3))
    res_all[[i]] = score_df[1:2,]
  }
  res_df = res_all[[1]]
  for(h in 2:length(res_all)){
    res_df = rbind(res_df, res_all[[h]])
  }
  res_df = res_df[!is.na(res_df$corr),]
  return(res_df)
}

ptm <- proc.time()


load("~/Dropbox/research/benchmark/rdata/9cellmix_qc.RData")
source('~/Dropbox/research/benchmark/trajectory_code/util_func.R')

sce_SC1_qc = filter_sce_genes(sce_SC1_qc)
sce_SC1_qc = scran_norm(sce_SC1_qc)
sce_SC1_qc = prep_traj_order(sce_SC1_qc)
sce_SC1_qc_traj = sce_SC1_qc[,colData(sce_SC1_qc)$traj =="YES"]
max_score_SC1 = slingshot_sup_order(sce_SC1_qc_traj,num_k=7)
max_score_SC1$method=method
max_score_SC1$design="cellmix"
max_score_SC1$dataset="cellmix1"

sce_SC2_qc = filter_sce_genes(sce_SC2_qc)
sce_SC2_qc = scran_norm(sce_SC2_qc)
sce_SC2_qc = prep_traj_order(sce_SC2_qc)
sce_SC2_qc_traj = sce_SC2_qc[,colData(sce_SC2_qc)$traj =="YES"]
max_score_SC2 = slingshot_sup_order(sce_SC2_qc_traj,num_k=7)
max_score_SC2$method=method
max_score_SC2$design="cellmix"
max_score_SC2$dataset="cellmix2"

sce_SC3_qc = filter_sce_genes(sce_SC3_qc)
sce_SC3_qc = scran_norm(sce_SC3_qc)
sce_SC3_qc = prep_traj_order(sce_SC3_qc)
sce_SC3_qc_traj = sce_SC3_qc[,colData(sce_SC3_qc)$traj =="YES"]
max_score_SC3 = slingshot_sup_order(sce_SC3_qc_traj,num_k=7)
max_score_SC3$method=method
max_score_SC3$design="cellmix"
max_score_SC3$dataset="cellmix3"

sce_SC4_qc = filter_sce_genes(sce_SC4_qc)
sce_SC4_qc = scran_norm(sce_SC4_qc)
sce_SC4_qc = prep_traj_order(sce_SC4_qc)
sce_SC4_qc_traj = sce_SC4_qc[,colData(sce_SC4_qc)$traj =="YES"]
max_score_SC4 = slingshot_sup_order(sce_SC4_qc_traj,num_k=7)
max_score_SC4$method=method
max_score_SC4$design="cellmix"
max_score_SC4$dataset="cellmix4"



load("~/Dropbox/research/benchmark/rdata/mRNAmix_qc.RData")

sce2_qc = filter_sce_genes(sce2_qc)
sce2_qc = scran_norm(sce2_qc)
sce2_qc = prep_RNA_traj_order(sce2_qc)
max_score_sce2 = slingshot_sup_order(sce2_qc,start_grp="1_0_0",num_k=7)
max_score_sce2$method=method
max_score_sce2$design="RNAmix"
max_score_sce2$dataset="RNAmix_CEL-seq2"

sce8_qc = filter_sce_genes(sce8_qc)
sce8_qc = scran_norm(sce8_qc)
sce8_qc = prep_RNA_traj_order(sce8_qc)
max_score_sce8 = slingshot_sup_order(sce8_qc,start_grp="1_0_0",num_k=7)
max_score_sce8$method=method
max_score_sce8$design="RNAmix"
max_score_sce8$dataset="RNAmix_Sort-seq"


final_res = rbind(max_score_SC1,max_score_SC2,max_score_SC3,max_score_SC4,
                  max_score_sce2, max_score_sce8)

write.csv(final_res,file=paste0("~/Dropbox/research/benchmark/traj_result/trajectory_result_",method,"_k_",7,"_.csv"),row.names = F)

final_time = (proc.time() - ptm)[3]

write.csv(data.frame(method=method,time_taken=final_time),file=paste0("~/Dropbox/research/benchmark/traj_result/trajectory_time_",method,"_.csv"),row.names = F)
