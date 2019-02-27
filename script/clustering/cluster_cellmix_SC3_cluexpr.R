
setwd("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit")
library(scran)
library(scater)
library(CellBench)
set_cellbench_threads(nthreads = 1)

log_file =paste("log_file/clustering_cellmix_SC3_clexpr",format(Sys.time(), "%a_%b_%d"),".txt",sep = ".")
cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "start preprocessing...\n"), file = log_file, append = TRUE)


# load all data 

# make a list of dataset for `Cellbench`


res2 = readRDS("rdata/final/cellmix_after_imputation.Rds")




scran_high_var = function(sce,topn=1000){
  var.fit <- trendVar(sce, method="loess", use.spikes=FALSE)
  var.out <- decomposeVar(sce, var.fit)
  hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
  rowData(sce)$hi_var = FALSE
  rowData(sce)$hi_var[rownames(rowData(sce)) %in% rownames(hvg.out)] = TRUE
  return(sce)
}




library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters="ensembl_gene_id",
                attributes=c("external_gene_name","ensembl_gene_id","external_gene_source"),
                values=rownames(res2$result[[1]]),
                mart=mart)

G_list = G_list[G_list$external_gene_source == "HGNC Symbol",]
G_list$format_name = paste("XXXX",G_list$external_gene_name,G_list$ensembl_gene_id,sep="_")
rownames(G_list) = G_list$ensembl_gene_id
format_rowname = function(sce){
  sce = sce[rownames(sce) %in% G_list$ensembl_gene_id,]
  rownames(sce) = G_list[rownames(sce),"format_name"]
  return(sce)
}

hivar_method = list(scran_hi = scran_high_var)
format_method = list(format_rowname=format_rowname)




# set the clustering methods


library(SC3)
library(clusterExperiment)
NUM_OF_THREAD=4


SC3_c = function(sce){
  tp = system.time({
    try_res = try({
      rowData(sce)$feature_symbol <- rownames(sce)
      sce = sce[rowData(sce)$hi_var,]
      sce <- sc3_estimate_k(sce)
      k_est <- sce@metadata$sc3$k_estimation
      sce <- sc3(sce, ks = k_est, biology = FALSE, n_cores=NUM_OF_THREAD,  k_estimator = FALSE, rand_seed=2333333)
      
      eval(parse(text=paste0("colData(sce)$clustering_res <- colData(sce)$sc3_", k_est, "_clusters")))
    })
    if (class(try_res) == "try-error") {
      cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
    }
  })
  
  method_name = "SC3"
  method_type = "clustering"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}


cluster_expr = function(sce){
  tp = system.time({
    try_res = try({
      sce = sce[rowData(sce)$hi_var,]
      se = RSEC(sce, isCount = FALSE, whichAssay = 2,minSizes=5,
                reduceMethod="PCA", nReducedDims=5,ncores=NUM_OF_THREAD, random.seed=176201)
      colData(sce)$clustering_res = primaryCluster(se)
    })
    if (class(try_res) == "try-error") {
      cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
    }
  })
  
  method_name = "clusterExperiment"
  method_type = "clustering"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}


clustering_method <- list(
  SC3 = SC3_c,
  clusterExperiment = cluster_expr
)



# apply clustering


cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "start to apply clustering methods...\n"), file = log_file, append = TRUE)

res2 = res2 %>% apply_methods(format_method)

res2 = res2 %>% apply_methods(hivar_method)

res2 = res2[,!(colnames(res2) %in% c("format_method","hivar_method"))]

res_c <- res2 %>%
  apply_methods(clustering_method)

cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "finish all clustering methods...\n"), file = log_file, append = TRUE)


clean_sce = function(sce){ # remove expr data after clustering to save space
  logcounts(sce) = NULL
  counts(sce) = NULL
  return(sce)
}
clean_up = list(cl = clean_sce)

res_c = res_c %>% apply_methods(clean_up)
res_c = res_c[,!(colnames(res_c) %in% c("clean_up"))]





saveRDS(res_c, file="rdata/cellmix_after_clustering_SC3_clexper.Rds")



# apply metrics for evaluation


library(mclust)


ARI_matric = function(sce){
  if(!("clustering_res" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    ari_val = adjustedRandIndex(sce$group, sce$clustering_res)
  }else{
    ari_val = adjustedRandIndex(sce$cell_line, sce$clustering_res)
  }
  
  return(ari_val)
}

clustering_evaluation <- list(
  ARI=ARI_matric
)

res_ARI <- res_c %>%
  apply_methods(clustering_evaluation)

cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "cellmix Done! saved the result to file.\n"), file = log_file, append = TRUE)



