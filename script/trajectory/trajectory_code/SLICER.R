library("SLICER")
library(lle)
method="SLICER"

get_max_score = function(col_anno, order_i, branches){
  
  tmp = table(branches[!is.na(col_anno$H2228_to_HCC827)])
  tmp = tmp[tmp>10]
  print(tmp)
  col_st = col_anno[branches %in% names(tmp),]
  cor1 = cor(order_i[branches %in% names(tmp)], col_st$H2228_to_HCC827, use = "pairwise.complete.obs")
  ov1 = sum(!is.na(col_st$H2228_to_HCC827))/sum(!is.na(col_anno$H2228_to_HCC827))
  
  tmp = table(branches[!is.na(col_anno$H1975_to_HCC827)])
  tmp = tmp[tmp>10]
  col_st = col_anno[branches %in% names(tmp),]
  cor2 = cor(order_i[branches %in% names(tmp)], col_st$H1975_to_HCC827, use = "pairwise.complete.obs")
  ov2 = sum(!is.na(col_st$H1975_to_HCC827))/sum(!is.na(col_anno$H1975_to_HCC827))
  
  tmp = table(branches[!is.na(col_anno$H2228_to_H1975)])
  tmp = tmp[tmp>10]
  col_st = col_anno[branches %in% names(tmp),]
  cor3 = cor(order_i[branches %in% names(tmp)], col_st$H2228_to_H1975, use = "pairwise.complete.obs")
  ov3 = sum(!is.na(col_st$H2228_to_H1975))/sum(!is.na(col_anno$H2228_to_H1975))
  
  res_df = data.frame(corr=abs(c(cor1,cor2,cor3)),overlap=c(ov1,ov2,ov3))
  res_df = res_df[order(res_df$overlap,decreasing = T),]
  print(res_df)
  return(res_df[1,])
}


slicer_order = function(sce, start_grp="9_0_0"){
  high_var_genes = scran_high_var(sce)
  res_all = list()
  t=1
  for(i in 1:10){
    set.seed(1000*i)
    pr = (1000:1)+100
    pr = pr/sum(pr)
    high_var_sel = sample(high_var_genes,size=500,prob=pr)
    dat = t(logcounts(sce)[high_var_sel,])
    res_df=tryCatch({
    k = select_k(dat, kmin = 2, kmax = 20, by = 3)
    traj_lle = lle(dat, m=2, k=k)$Y
    traj_graph = conn_knn_graph(traj_lle,5)
    extreme_cells = find_extreme_cells(traj_graph, traj_lle)
    if("9_0_0" %in% colData(sce)$group[extreme_cells]){
      H2228_st = extreme_cells[(colData(sce)$group == start_grp)[extreme_cells]][1]
      cells_ordered = cell_order(traj_graph, H2228_st)
      branches = assign_branches(traj_graph,H2228_st)
      res_df = get_max_score(as.data.frame(colData(sce)), cells_ordered)
    }else{
      other_st = extreme_cells[1]
      cells_ordered = cell_order(traj_graph, other_st)
      branches = assign_branches(traj_graph,other_st)
      res_df = get_max_score(as.data.frame(colData(sce)), cells_ordered, branches)
    }
    res_df}, error=function(cond){print(cond);return(NA)})
    if(!is.na(res_df)){
      res_df[order(res_df$corr,decreasing = T),]
      res_df = res_df[1,]
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
max_score_SC1 = slicer_order(sce_SC1_qc_traj)
max_score_SC1$method=method
max_score_SC1$design="cellmix"
max_score_SC1$dataset="cellmix1"


sce_SC2_qc = filter_sce_genes(sce_SC2_qc)
sce_SC2_qc = scran_norm(sce_SC2_qc)
sce_SC2_qc = prep_traj_order(sce_SC2_qc)
sce_SC2_qc_traj = sce_SC2_qc[,colData(sce_SC2_qc)$traj =="YES"]
max_score_SC2 = slicer_order(sce_SC2_qc_traj)
max_score_SC2$method=method
max_score_SC2$design="cellmix"
max_score_SC2$dataset="cellmix2"

sce_SC3_qc = filter_sce_genes(sce_SC3_qc)
sce_SC3_qc = scran_norm(sce_SC3_qc)
sce_SC3_qc = prep_traj_order(sce_SC3_qc)
sce_SC3_qc_traj = sce_SC3_qc[,colData(sce_SC3_qc)$traj =="YES"]
max_score_SC3 = slicer_order(sce_SC3_qc_traj)
max_score_SC3$method=method
max_score_SC3$design="cellmix"
max_score_SC3$dataset="cellmix3"

sce_SC4_qc = filter_sce_genes(sce_SC4_qc)
sce_SC4_qc = scran_norm(sce_SC4_qc)
sce_SC4_qc = prep_traj_order(sce_SC4_qc)
sce_SC4_qc_traj = sce_SC4_qc[,colData(sce_SC4_qc)$traj =="YES"]
max_score_SC4 = slicer_order(sce_SC4_qc_traj)
max_score_SC4$method=method
max_score_SC4$design="cellmix"
max_score_SC4$dataset="cellmix4"



load("~/Dropbox/research/benchmark/rdata/mRNAmix_qc.RData")

sce2_qc = filter_sce_genes(sce2_qc)
sce2_qc = scran_norm(sce2_qc)
sce2_qc = prep_RNA_traj_order(sce2_qc)
max_score_sce2 = slicer_order(sce2_qc,start_grp="1_0_0")
max_score_sce2$method=method
max_score_sce2$design="RNAmix"
max_score_sce2$dataset="RNAmix_CEL-seq2"

sce8_qc = filter_sce_genes(sce8_qc)
sce8_qc = scran_norm(sce8_qc)
sce8_qc = prep_RNA_traj_order(sce8_qc)
max_score_sce8 = slicer_order(sce8_qc,start_grp="1_0_0")
max_score_sce8$method=method
max_score_sce8$design="RNAmix"
max_score_sce8$dataset="RNAmix_Sort-seq"


final_res = rbind(max_score_SC1,max_score_SC2,max_score_SC3,max_score_SC4,
                  max_score_sce2, max_score_sce8)

write.csv(final_res,file=paste0("~/Dropbox/research/benchmark/traj_result/trajectory_result_",method,"_.csv"),row.names = F)

final_time = (proc.time() - ptm)[3]

write.csv(data.frame(method=method,time_taken=final_time),file=paste0("~/Dropbox/research/benchmark/traj_result/trajectory_time_",method,"_.csv"),row.names = F)

# 
# high_var_genes = scran_high_var(sce_SC1_qc_traj)
# #genes = select_genes(traj)
# k = select_k(t(logcounts(sce_SC1_qc_traj[high_var_genes,])), kmin=2)
# traj_lle = lle(t(logcounts(sce_SC1_qc_traj[high_var_genes,])), m=2, k=k)$Y
# traj_graph = conn_knn_graph(traj_lle,5)
# ends = find_extreme_cells(traj_graph, traj_lle)
# start = 142
# cells_ordered = cell_order(traj_graph, start)
# branches = assign_branches(traj_graph,start)
# 
# 
# ggplot(data=NULL,aes(x=traj_lle[,1],y=traj_lle[,2],col=factor(branches)))+geom_point()
