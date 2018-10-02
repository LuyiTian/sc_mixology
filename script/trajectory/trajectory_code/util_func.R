library(scran)
library(scater)

prep_traj_order = function(sce){
  sce$group = paste(colData(sce)$H2228,colData(sce)$H1975,colData(sce)$HCC827,sep="_")
  sce$H2228_to_H1975 = NA
  sce$H2228_to_H1975[sce$group=="9_0_0"] = 0
  sce$H2228_to_H1975[sce$group=="7_1_1"] = 1
  sce$H2228_to_H1975[sce$group=="5_2_2"] = 2
  sce$H2228_to_H1975[sce$group=="3_3_3"] = 3
  sce$H2228_to_H1975[sce$group=="2_5_2"] = 4
  sce$H2228_to_H1975[sce$group=="1_7_1"] = 5
  sce$H2228_to_H1975[sce$group=="0_9_0"] = 6
  
  sce$H2228_to_HCC827 = NA
  sce$H2228_to_HCC827[sce$group=="9_0_0"] = 0
  sce$H2228_to_HCC827[sce$group=="7_1_1"] = 1
  sce$H2228_to_HCC827[sce$group=="5_2_2"] = 2
  sce$H2228_to_HCC827[sce$group=="3_3_3"] = 3
  sce$H2228_to_HCC827[sce$group=="2_2_5"] = 4
  sce$H2228_to_HCC827[sce$group=="1_1_7"] = 5
  sce$H2228_to_HCC827[sce$group=="0_0_9"] = 6
  
  sce$H1975_to_HCC827 = NA
  sce$H1975_to_HCC827[sce$group=="0_9_0"] = 0
  sce$H1975_to_HCC827[sce$group=="1_7_1"] = 1
  sce$H1975_to_HCC827[sce$group=="2_5_2"] = 2
  sce$H1975_to_HCC827[sce$group=="3_3_3"] = 3
  sce$H1975_to_HCC827[sce$group=="2_2_5"] = 4
  sce$H1975_to_HCC827[sce$group=="1_1_7"] = 5
  sce$H1975_to_HCC827[sce$group=="0_0_9"] = 6
  return(sce)
}



prep_RNA_traj_order = function(sce){
  sce$group = paste(colData(sce)$H2228_prop,colData(sce)$H1975_prop,colData(sce)$HCC827_prop,sep="_")
  sce$H2228_to_H1975 = NA
  sce$H2228_to_H1975[sce$group=="1_0_0"] = 0
  sce$H2228_to_H1975[sce$group=="0.68_0.16_0.16"] = 1
  sce$H2228_to_H1975[sce$group=="0.33_0.33_0.33"] = 2
  sce$H2228_to_H1975[sce$group=="0.16_0.68_0.16"] = 3
  sce$H2228_to_H1975[sce$group=="0_1_0"] = 4
  
  sce$H2228_to_HCC827 = NA
  sce$H2228_to_HCC827[sce$group=="1_0_0"] = 0
  sce$H2228_to_HCC827[sce$group=="0.68_0.16_0.16"] = 1
  sce$H2228_to_HCC827[sce$group=="0.33_0.33_0.33"] = 2
  sce$H2228_to_HCC827[sce$group=="0.16_0.16_0.68"] = 3
  sce$H2228_to_HCC827[sce$group=="0_0_1"] = 4
  
  sce$H1975_to_HCC827 = NA
  sce$H1975_to_HCC827[sce$group=="0_1_0"] = 0
  sce$H1975_to_HCC827[sce$group=="0.16_0.68_0.16"] = 1
  sce$H1975_to_HCC827[sce$group=="0.33_0.33_0.33"] = 2
  sce$H1975_to_HCC827[sce$group=="0.16_0.16_0.68"] = 3
  sce$H1975_to_HCC827[sce$group=="0_0_1"] = 4
  return(sce)
}



filter_sce_genes = function(sce){
  keep1 = (apply(counts(sce), 1, function(x) mean(x[x>0])) > 1)  # average count larger than one
  keep2 = (rowSums(counts(sce)>0) > 0.05*ncol(sce))  # expressed in at least six cells
  sce = sce[(keep1 & keep2), ]
  return(sce)
}

scran_norm = function(sce){
  sce = computeSumFactors(sce)
  sce = normalize(sce)
  return(sce)
}


scran_high_var = function(sce,topn=1000){
  var.fit <- trendVar(sce, method="loess", use.spikes=FALSE)
  var.out <- decomposeVar(sce, var.fit)
  hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
  return(rownames(hvg.out))
}
