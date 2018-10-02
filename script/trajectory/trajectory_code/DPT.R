library(dpt)
library(destiny)
method="DPT"
library(scran)
library(scater)

get_max_score = function(col_anno, pseudotime, branch){
  
  cor1 = cor(pseudotime, col_anno$H2228_to_H1975, use = "pairwise.complete.obs")
  cor2 = cor(pseudotime, col_anno$H2228_to_HCC827, use = "pairwise.complete.obs")
  
  ov1 = max(sum(branch[!is.na(col_anno$H2228_to_H1975)] %in% c("branch 1","branch 2"))/sum(!is.na(col_anno$H2228_to_H1975)),
            sum(branch[!is.na(col_anno$H2228_to_H1975)] %in% c("branch 1","branch 3"))/sum(!is.na(col_anno$H2228_to_H1975)),
            sum(branch[!is.na(col_anno$H2228_to_H1975)] %in% c("branch 2","branch 3"))/sum(!is.na(col_anno$H2228_to_H1975)))
  ov2 = max(sum(branch[!is.na(col_anno$H2228_to_HCC827)] %in% c("branch 1","branch 2"))/sum(!is.na(col_anno$H2228_to_HCC827)),
            sum(branch[!is.na(col_anno$H2228_to_HCC827)] %in% c("branch 1","branch 3"))/sum(!is.na(col_anno$H2228_to_HCC827)),
            sum(branch[!is.na(col_anno$H2228_to_HCC827)] %in% c("branch 2","branch 3"))/sum(!is.na(col_anno$H2228_to_HCC827)))
  
  res_df = data.frame(corr=abs(c(cor1,cor2)),overlap=c(ov1,ov2))
  #print(res_df)
  res_df = res_df[order(res_df$overlap,decreasing = T),]
  return(res_df)
}



DPT_order = function(sce, start_grp="9_0_0"){
  res_all=list()
  high_var_genes = scran_high_var(sce)
  for(i in 1:10){
    set.seed(1000*i) # so the result is reproducible
    pr = (1000:1)+100
    pr = pr/sum(pr)
    high_var_sel = sample(high_var_genes,size=500,prob=pr)
    sig = find.sigmas(t(logcounts(sce)[high_var_sel,]))
    ts <- Transitions(t(logcounts(sce)[high_var_sel,]), sig,k=10)
    pt <- dpt(ts, branching = TRUE,root=which(sce$group == start_grp)[1])
    score_df = get_max_score(as.data.frame(colData(sce)),pt$DPT, pt$Branch)
    res_all[[i]] = score_df
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
max_score_SC1 = DPT_order(sce_SC1_qc_traj)
max_score_SC1$method=method
max_score_SC1$design="cellmix"
max_score_SC1$dataset="cellmix1"


sce_SC2_qc = filter_sce_genes(sce_SC2_qc)
sce_SC2_qc = scran_norm(sce_SC2_qc)
sce_SC2_qc = prep_traj_order(sce_SC2_qc)
sce_SC2_qc_traj = sce_SC2_qc[,colData(sce_SC2_qc)$traj =="YES"]
max_score_SC2 = DPT_order(sce_SC2_qc_traj)
max_score_SC2$method=method
max_score_SC2$design="cellmix"
max_score_SC2$dataset="cellmix2"

sce_SC3_qc = filter_sce_genes(sce_SC3_qc)
sce_SC3_qc = scran_norm(sce_SC3_qc)
sce_SC3_qc = prep_traj_order(sce_SC3_qc)
sce_SC3_qc_traj = sce_SC3_qc[,colData(sce_SC3_qc)$traj =="YES"]
max_score_SC3 = DPT_order(sce_SC3_qc_traj)
max_score_SC3$method=method
max_score_SC3$design="cellmix"
max_score_SC3$dataset="cellmix3"

sce_SC4_qc = filter_sce_genes(sce_SC4_qc)
sce_SC4_qc = scran_norm(sce_SC4_qc)
sce_SC4_qc = prep_traj_order(sce_SC4_qc)
sce_SC4_qc_traj = sce_SC4_qc[,colData(sce_SC4_qc)$traj =="YES"]
max_score_SC4 = DPT_order(sce_SC4_qc_traj)
max_score_SC4$method=method
max_score_SC4$design="cellmix"
max_score_SC4$dataset="cellmix4"



load("~/Dropbox/research/benchmark/rdata/mRNAmix_qc.RData")

sce2_qc = filter_sce_genes(sce2_qc)
sce2_qc = scran_norm(sce2_qc)
sce2_qc = prep_RNA_traj_order(sce2_qc)
max_score_sce2 = DPT_order(sce2_qc,start_grp="1_0_0")
max_score_sce2$method=method
max_score_sce2$design="RNAmix"
max_score_sce2$dataset="RNAmix_CEL-seq2"

sce8_qc = filter_sce_genes(sce8_qc)
sce8_qc = scran_norm(sce8_qc)
sce8_qc = prep_RNA_traj_order(sce8_qc)
max_score_sce8 = DPT_order(sce8_qc,start_grp="1_0_0")
max_score_sce8$method=method
max_score_sce8$design="RNAmix"
max_score_sce8$dataset="RNAmix_Sort-seq"


final_res = rbind(max_score_SC1,max_score_SC2,max_score_SC3,max_score_SC4,
                  max_score_sce2, max_score_sce8)

write.csv(final_res,file=paste0("~/Dropbox/research/benchmark/traj_result/trajectory_result_",method,"_.csv"),row.names = F)

final_time = (proc.time() - ptm)[3]

write.csv(data.frame(method=method,time_taken=final_time),file=paste0("~/Dropbox/research/benchmark/traj_result/trajectory_time_",method,"_.csv"),row.names = F)

