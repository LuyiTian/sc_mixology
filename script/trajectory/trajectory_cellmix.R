

setwd("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit")
library(scran)
library(scater)
library(CellBench)
library(R.utils)
set_cellbench_threads(nthreads = 1)
MAX_TIME = 60*120
NUM_OF_THREAD=8

log_file =paste("log_file/trajectory_cellmix",format(Sys.time(), "%a_%b_%d"),"txt",sep = ".")
cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "start preprocessing...\n"), file = log_file, append = TRUE)



# load the data
res2 = readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/cellmix_all_after_imputation.Rds")


# feature selection
scran_high_var = function(sce,topn=1000){
  if(max(logcounts(sce)>100)){
    logcounts(sce) = log2(logcounts(sce)+1)
  }
  var.fit <- trendVar(sce, method="loess", use.spikes=FALSE)
  var.out <- decomposeVar(sce, var.fit)
  hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
  rowData(sce)$hi_var = FALSE
  rowData(sce)$hi_var[rownames(rowData(sce)) %in% rownames(hvg.out)] = TRUE
  return(sce)
}
hivar_method = list(scran_hi = scran_high_var)

# get trajectory annotation
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

prep_traj_warpper = function(sce){
  if("9 0 0" %in% sce$group){
    return(prep_traj_order(sce))
  }else{
    return(prep_RNA_traj_order(sce))
  }
}

prep_traj <- list(
  prep_traj=prep_traj_warpper
)


#####
##### trajectory method

library(monocle)
library(TSCAN)
library(mclust, quietly = TRUE)
library(destiny)
library(dpt)
library(slingshot)
library(SLICER)
library(lle)

monocle2_DDRTree = function(sce){
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
  if("9 0 0" %in% sce$group){
    start_grp="9 0 0"
  }else{
    start_grp="1 0 0"
  }
  cat(paste(format(Sys.time(), "%a %b %d %X %Y. START monocle2_DDRTree: "), print(paste(metadata(sce)$running_time$method,sep=",",collapse=",")),"\n"), file = log_file, append = TRUE)
  tp = system.time({
    try_res = try({
      withTimeout({
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
      }
      tmp = unlist(lapply(max_score_df,function(x){colMeans(x)[2]}))
      metadata(sce)$traj_eval = max_score_df[[which(tmp==max(tmp))]]
      },timeout = MAX_TIME, onTimeout="error")
    })
  })
  if (class(try_res) == "try-error") {
    cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
  }
  
  
  method_name="monocle2_DDRTree"
  method_type="trajectory"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  
  return(sce)
  
}



monocle2_SimplePPT = function(sce){
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
  if("9 0 0" %in% sce$group){
    start_grp="9 0 0"
  }else{
    start_grp="1 0 0"
  }
  cat(paste(format(Sys.time(), "%a %b %d %X %Y. START monocle2_SimplePPT: "), print(paste(metadata(sce)$running_time$method,sep=",",collapse=",")),"\n"), file = log_file, append = TRUE)
  tp = system.time({
    try_res = try({
      withTimeout({
        pd = new("AnnotatedDataFrame", data = as.data.frame(colData(sce)))
        cds = new("CellDataSet", exprs = logcounts(sce), phenoData = pd, expressionFamily = VGAM::negbinomial.size())
        sizeFactors(cds)=1 # assume the data has been normalized.
        if("hi_var" %in% colnames(rowData(sce))){
          cds = setOrderingFilter(cds, ordering_genes = rownames(sce)[rowData(sce)$hi_var])
        }
        cds = reduceDimension(cds, method = "SimplePPT",norm_method="none",pseudo_expr=0)
        cds = orderCells(cds)
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
        }
        tmp = unlist(lapply(max_score_df,function(x){colMeans(x)[2]}))
        metadata(sce)$traj_eval = max_score_df[[which(tmp==max(tmp))]]
      },timeout = MAX_TIME, onTimeout="error")
    })
  })
  if (class(try_res) == "try-error") {
    cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
  }
  
  
  method_name="monocle2_SimplePPT"
  method_type="trajectory"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  
  return(sce)
  
}



TSCAN_order = function(sce){
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
  if("9 0 0" %in% sce$group){
    start_grp="9 0 0"
  }else{
    start_grp="1 0 0"
  }
  cat(paste(format(Sys.time(), "%a %b %d %X %Y. START TSCAN: "), print(paste(metadata(sce)$running_time$method,sep=",",collapse=",")),"\n"), file = log_file, append = TRUE)
  tp = system.time({
    try_res = try({
      withTimeout({
      lpsmclust = exprmclust(logcounts(sce)[rowData(sce)$hi_var,],
                             clusternum = 10)
      tmp = table(lpsmclust$clusterid[colData(sce)$group==start_grp])  # specify H2228 as root state.
      tmp = tmp[order(tmp,decreasing = T)]
      lpsorder <- TSCANorder(lpsmclust,orderonly=F,listbranch=T)
      corr = c()
      ov = c()
      for(j in 1:length(lpsorder)){
        if(lpsorder[[j]]$State[1] %in% names(tmp[1]) | lpsorder[[j]]$State[length(lpsorder[[j]]$State)] %in% names(tmp[1])){
          max_j = get_max_score(as.data.frame(colData(sce)), lpsorder[[j]])
          corr = c(corr, max_j$corr)
          ov = c(ov, max_j$overlap)
        }
      }
      res_df = data.frame(corr=corr,overlap=ov)
      res_df[order(res_df$corr,decreasing = T),]
      res_df = res_df[1:2,]
      metadata(sce)$traj_eval = res_df
    },timeout = MAX_TIME, onTimeout="error")
    })
  })
  if (class(try_res) == "try-error") {
    cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
  }
  
  method_name = "TSCAN"
  method_type = "trajectory"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}



slingshot_order = function(sce){
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
  if("9 0 0" %in% sce$group){
    start_grp="9 0 0"
  }else{
    start_grp="1 0 0"
  }
  cat(paste(format(Sys.time(), "%a %b %d %X %Y. START slingshot: "), print(paste(metadata(sce)$running_time$method,sep=",",collapse=",")),"\n"), file = log_file, append = TRUE)
  tp = system.time({
    try_res = try({
      withTimeout({
      sce_de_traj = plotPCA(sce,return_SCE=T, draw_plot=T, rerun=T,run_args=list(feature_set=rownames(sce)[rowData(sce)$hi_var]))
      #kmeans_clu = kmeans(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"),centers=num_k,iter.max = 10000)
      cl1 <- Mclust(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"))$classification
      colData(sce_de_traj)$GMM_cluster= as.factor(cl1)
      tmp = table(colData(sce_de_traj)$GMM_cluster[colData(sce_de_traj)$group==start_grp])  # specify H2228 as root state.
      tmp = tmp[order(tmp,decreasing = T)]
      
      slingshot_lin <- slingshot::getLineages(SingleCellExperiment::reducedDim(sce_de_traj,"PCA"), colData(sce_de_traj)$GMM_cluster,start.clus = names(tmp)[1])
      slingshot_crv <- slingshot::getCurves(slingshot_lin)
      slingshot_pseudo <- slingshot::pseudotime(slingshot_crv)
      col_data_df = as.data.frame(colData(sce_de_traj))
      score_df = apply(slingshot_pseudo,2,function(x){get_max_score(col_data_df, x)})
      score_df = Reduce(rbind,score_df)
      score_df = score_df[order(score_df$corr,decreasing = T),]
      #print(lines(slingshot_crv, lwd = 3))
      metadata(sce)$traj_eval = score_df[1:2,]
      },timeout = MAX_TIME, onTimeout="error")
    })
  })
  if (class(try_res) == "try-error") {
    cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
  }
  
  method_name = "Slingshot"
  method_type = "trajectory"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}


DPT_order = function(sce){
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
  if("9 0 0" %in% sce$group){
    start_grp="9 0 0"
  }else{
    start_grp="1 0 0"
  }
  cat(paste(format(Sys.time(), "%a %b %d %X %Y. START DPT: "), print(paste(metadata(sce)$running_time$method,sep=",",collapse=",")),"\n"), file = log_file, append = TRUE)
  tp = system.time({
    try_res = try({
      withTimeout({
        sig = find.sigmas(t(logcounts(sce)[rowData(sce)$hi_var,]))
        ts <- Transitions(t(logcounts(sce)[rowData(sce)$hi_var,]), sig,k=10)
        pt <- dpt(ts, branching = TRUE,root=which(sce$group == start_grp)[1])
        score_df = get_max_score(as.data.frame(colData(sce)),pt$DPT, pt$Branch)
        metadata(sce)$traj_eval = score_df
    },timeout = MAX_TIME, onTimeout="error")
    })
  })
  if (class(try_res) == "try-error") {
    cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
  }
  
  method_name = "DPT"
  method_type = "trajectory"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}



slicer_order = function(sce){
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
  if("9 0 0" %in% sce$group){
    start_grp="9 0 0"
  }else{
    start_grp="1 0 0"
  }
  cat(paste(format(Sys.time(), "%a %b %d %X %Y. START slicer: "), print(paste(colnames(metadata(sce)$running_time),sep=",",collapse=",")),"\n"), file = log_file, append = TRUE)
  tp = system.time({
   try_res = try({
      withTimeout({
        dat = t(logcounts(sce)[rowData(sce)$hi_var,])
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
        res_df[order(res_df$corr,decreasing = T),]
        res_df = res_df[1,]
        metadata(sce)$traj_eval = res_df
    },timeout = MAX_TIME, onTimeout="error")
   })
  })
  if (class(try_res) == "try-error") {
    cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
  }
  
  method_name = "SLICER"
  method_type = "trajectory"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}




get_trajectory_corr = function(sce){
  if("traj_eval" %in% metadata(sce)){
    return(metadata(sce)$traj_eval$corr)
  }else{
    return(NA)
  }
}

get_trajectory_overlap = function(sce){
  if("traj_eval" %in% metadata(sce)){
    return(metadata(sce)$traj_eval$overlap)
  }else{
    return(NA)
  }
}



traj_method <- list(
  Monocle2_DDRTree=monocle2_DDRTree,
  TSCAN=TSCAN_order,
  DPT=DPT_order,
  Slingshot=slingshot_order,
  SLICER=slicer_order
)





cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "start to apply trajectory methods...\n"), file = log_file, append = TRUE)
res2 = res2[unlist(lapply(res2$result, function(x){"logcounts" %in% assayNames(x)})),]
res2 = res2[!(res2$impute_method=="knn_smooth2"),] # trajectory methods would clash on knn_smooth2 generated matrix.
res2 = res2 %>% apply_methods(hivar_method)

res2 = res2 %>% apply_methods(prep_traj)

res2 = res2[,!(colnames(res2) %in% c("prep_traj","hivar_method"))]


res_c <- res2 %>%
  apply_methods(traj_method)

cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "finish all trajectory methods...\n"), file = log_file, append = TRUE)


### post clean up to save space

clean_sce = function(sce){ # remove expr data after clustering to save space
  logcounts(sce) = NULL
  counts(sce) = NULL
  return(sce)
}
clean_up = list(cl = clean_sce)

res_c = res_c %>% apply_methods(clean_up)
res_c = res_c[,!(colnames(res_c) %in% c("clean_up"))]



saveRDS(res_c, file="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/cellmix_all_after_trajectory.Rds")



