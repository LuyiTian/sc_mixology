# TSCAN
library(TSCAN)
method="TSCAN"
get_max_score = function(col_anno, order_i){
  col_st = col_anno[order_i$sample_name,]
  
  cor1 = cor(order_i$Pseudotime, col_st$H2228_to_H1975, use = "pairwise.complete.obs")
  cor2 = cor(order_i$Pseudotime, col_st$H2228_to_HCC827, use = "pairwise.complete.obs")
  
  ov1 = sum(!is.na(col_st$H2228_to_H1975))/sum(!is.na(col_anno$H2228_to_H1975))
  ov2 = sum(!is.na(col_st$H2228_to_HCC827))/sum(!is.na(col_anno$H2228_to_HCC827))
  
  res_df = data.frame(corr=abs(c(cor1,cor2)),overlap=c(ov1,ov2))
  res_df = res_df[order(res_df$overlap,decreasing = T),]
  #print(res_df)
  return(res_df[1,])
}

tscan_order = function(sce, start_grp="9_0_0"){
  high_var_genes = scran_high_var(sce)
  res_all = list()
  t=1
  for(i in 1:10){
    set.seed(1000*i)
    pr = (1000:1)+100
    pr = pr/sum(pr)
    high_var_sel = sample(high_var_genes,size=500,prob=pr)
    res_df=tryCatch({
    lpsmclust = exprmclust(logcounts(sce)[high_var_sel,],
                           clusternum = 10)
    tmp = table(lpsmclust$clusterid[colData(sce)$group==start_grp])  # specify H2228 as root state.
    tmp = tmp[order(tmp,decreasing = T)]
    lpsorder <- TSCANorder(lpsmclust,orderonly=F,listbranch=T)
    corr = c()
    ov = c()
    for(j in 1:length(lpsorder)){
      #print(unique(lpsorder[[j]]$State))
      if(lpsorder[[j]]$State[1] %in% names(tmp[1]) | lpsorder[[j]]$State[length(lpsorder[[j]]$State)] %in% names(tmp[1])){
        max_j = get_max_score(as.data.frame(colData(sce)), lpsorder[[j]])
        corr = c(corr, max_j$corr)
        ov = c(ov, max_j$overlap)
      }
    }
    res_df = data.frame(corr=corr,overlap=ov)
    #print(res_df)
    res_df[order(res_df$corr,decreasing = T),]
    res_df = res_df[1:2,]
    }, error=function(cond){return(NA)}
    )
    if(!is.na(res_df)){
      res_all[[t]] = res_df
      t = t+1
    }
  }
  res_df = res_all[[1]]
  for(h in 2:length(res_all)){
    res_df = rbind(res_df, res_all[[h]])
  }
  return(res_df)
}


ptm <- proc.time()

load("~/Dropbox/research/benchmark/rdata/9cellmix_qc.RData")
source('~/Dropbox/research/benchmark/trajectory_code/util_func.R')

sce_SC1_qc = filter_sce_genes(sce_SC1_qc)
sce_SC1_qc = scran_norm(sce_SC1_qc)
sce_SC1_qc = prep_traj_order(sce_SC1_qc)
sce_SC1_qc_traj = sce_SC1_qc[,colData(sce_SC1_qc)$traj =="YES"]
max_score_SC1 = tscan_order(sce_SC1_qc_traj)
max_score_SC1$method=method
max_score_SC1$design="cellmix"
max_score_SC1$dataset="cellmix1"


sce_SC2_qc = filter_sce_genes(sce_SC2_qc)
sce_SC2_qc = scran_norm(sce_SC2_qc)
sce_SC2_qc = prep_traj_order(sce_SC2_qc)
sce_SC2_qc_traj = sce_SC2_qc[,colData(sce_SC2_qc)$traj =="YES"]
max_score_SC2 = tscan_order(sce_SC2_qc_traj)
max_score_SC2$method=method
max_score_SC2$design="cellmix"
max_score_SC2$dataset="cellmix1"

sce_SC3_qc = filter_sce_genes(sce_SC3_qc)
sce_SC3_qc = scran_norm(sce_SC3_qc)
sce_SC3_qc = prep_traj_order(sce_SC3_qc)
sce_SC3_qc_traj = sce_SC3_qc[,colData(sce_SC3_qc)$traj =="YES"]
max_score_SC3 = tscan_order(sce_SC3_qc_traj)
max_score_SC3$method=method
max_score_SC3$design="cellmix"
max_score_SC3$dataset="cellmix1"

sce_SC4_qc = filter_sce_genes(sce_SC4_qc)
sce_SC4_qc = scran_norm(sce_SC4_qc)
sce_SC4_qc = prep_traj_order(sce_SC4_qc)
sce_SC4_qc_traj = sce_SC4_qc[,colData(sce_SC4_qc)$traj =="YES"]
max_score_SC4 = tscan_order(sce_SC4_qc_traj)
max_score_SC4$method=method
max_score_SC4$design="cellmix"
max_score_SC4$dataset="cellmix1"



load("~/Dropbox/research/benchmark/rdata/mRNAmix_qc.RData")

sce2_qc = filter_sce_genes(sce2_qc)
sce2_qc = scran_norm(sce2_qc)
sce2_qc = prep_RNA_traj_order(sce2_qc)
max_score_sce2 = tscan_order(sce2_qc,start_grp="1_0_0")
max_score_sce2$method=method
max_score_sce2$design="RNAmix"
max_score_sce2$dataset="RNAmix_CEL-seq2"

sce8_qc = filter_sce_genes(sce8_qc)
sce8_qc = scran_norm(sce8_qc)
sce8_qc = prep_RNA_traj_order(sce8_qc)
max_score_sce8 = tscan_order(sce8_qc,start_grp="1_0_0")
max_score_sce8$method=method
max_score_sce8$design="RNAmix"
max_score_sce8$dataset="RNAmix_Sort-seq"




final_res = rbind(max_score_SC1,max_score_SC2,max_score_SC3,max_score_SC4,
                  max_score_sce2, max_score_sce8)

write.csv(final_res,file=paste0("~/Dropbox/research/benchmark/traj_result/trajectory_result_",method,".csv"),row.names = F)



final_time = (proc.time() - ptm)[3]

write.csv(data.frame(method=method,time_taken=final_time),file=paste0("~/Dropbox/research/benchmark/traj_result/trajectory_time_",method,"_.csv"),row.names = F)




