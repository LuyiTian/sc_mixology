# norm evaluation
setwd("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit")

library(CellBench)
library(ggplot2)
library(scater)
library(scran)
library(tidyr)
library(plyr)
library(cluster)
library(ggpubr)

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

silhouette_pca = function(sce){
  if("logcounts" %in% assayNames(sce)){
    try_res = try({
      sce = runPCA(sce,exprs_values="logcounts")
    })
    if (class(try_res) == "try-error") {
      return(NA)
    }
  }else{
    return(NA)
  }
  ret_val = NA
  try_res = try({
    if("cell_line_demuxlet" %in% colnames(colData(sce))){
      sil = silhouette(as.numeric(factor(sce$cell_line_demuxlet)),  dist(reducedDim(sce,"PCA")))
      ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))
    }else if("group" %in% colnames(colData(sce))){
      sil = silhouette(as.numeric(factor(sce$group)),  dist(reducedDim(sce,"PCA")))
      ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))    
    }else{
      ret_val=NA
    }
  })
  if (class(try_res) == "try-error") {
    return(NA)
  }
  return(ret_val)
}

silhouette_raw = function(sce){
  if("logcounts" %in% assayNames(sce)){
    try_res = try({
      sce = runPCA(sce,exprs_values="logcounts")
    })
    if (class(try_res) == "try-error") {
      return(NA)
    }
  }else{
    return(NA)
  }
  ret_val = NA
  try_res = try({
    if("cell_line_demuxlet" %in% colnames(colData(sce))){
      sil = silhouette(as.numeric(factor(sce$cell_line_demuxlet)),  dist(t(logcounts(sce)[rowData(sce)$hi_var,])))
      ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))
    }else if("group" %in% colnames(colData(sce))){
      sil = silhouette(as.numeric(factor(sce$group)),    dist(t(logcounts(sce)[rowData(sce)$hi_var,])))
      ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))    
    }else{
      ret_val=NA
    }
  })
  if (class(try_res) == "try-error") {
    return(NA)
  }
  return(ret_val)
}

#Association of expression measures with factors of wanted and unwanted variation.

fac_corr_u = function(sce){
  if("logcounts" %in% assayNames(sce)){
    try_res = try({
      sce = runPCA(sce,exprs_values="logcounts",ncomponents = 5)
    })
    if (class(try_res) == "try-error") {
      return(NA)
    }
  }else{
    return(NA)
  }
  score = 0
  pca_res = reducedDim(sce)
  vec = log10(sce$total_count_per_cell)
  for (i in 1:dim(pca_res)[2]){
    r = summary(lm(pca_res[,i]~vec))$r.squared
    score = c(score,r*attr(pca_res,"percentVar")[i])
  }
  ret_val = sum(score)
  return(ret_val)
}

fac_corr_w = function(sce){
  if("logcounts" %in% assayNames(sce)){
    try_res = try({
      sce = runPCA(sce,exprs_values="logcounts",ncomponents = 5)
    })
    if (class(try_res) == "try-error") {
      return(NA)
    }
  }else{
    return(NA)
  }
  score = 0
  pca_res = reducedDim(sce)
  if("group" %in% colnames(colData(sce))){
    vec = sce$group
  }else{
    vec = sce$cell_line_demuxlet
  }
  for (i in 1:dim(pca_res)[2]){
    r = summary(lm(pca_res[,i]~vec))$r.squared
    score = c(score,r*attr(pca_res,"percentVar")[i])
  }
  ret_val = sum(score)
  return(ret_val)
}

cor_RNAmix_within <- function(sce) {
  cor_val = c()
  norm_mat = logcounts(sce)
  if (max(norm_mat)>100){
    norm_mat = log2(norm_mat+1-min(0,min(norm_mat)))
  }
  for(group in unique(sce$group)){
    cor_mat = cor(norm_mat[,sce$group==group])
    cor_mat[!lower.tri(cor_mat)] = NA
    cor_mat = as.numeric(cor_mat)
    cor_mat = cor_mat[!is.na(cor_mat)]
    cor_val = c(cor_val,cor_mat)
  }
  return(mean(cor_val))
}

get_cor_two_groups = function(sce,grp1,grp2){
  norm_mat = logcounts(sce)
  if (max(norm_mat)>100){
    norm_mat = log2(norm_mat+1-min(0,min(norm_mat)))
  }
  da1 = norm_mat[,sce$group %in% grp1]
  da2 = norm_mat[,sce$group %in% grp2]
  cor_mat = cor(da1,da2)
  cor_mat[!lower.tri(cor_mat)] = NA
  cor_mat = as.numeric(cor_mat)
  cor_mat = cor_mat[!is.na(cor_mat)]
  return(cor_mat)
}

cor_RNAmix_between <- function(sce) {
  cor_val = c()
  group="0 0 1"
  grp2=c("0.16 0.16 0.68")
  cor_val = c(cor_val,get_cor_two_groups(sce,group,grp2))
  
  group="0 1 0"
  grp2=c("0.16 0.68 0.16")
  cor_val = c(cor_val,get_cor_two_groups(sce,group,grp2))
  
  group="1 0 0"
  grp2="0.68 0.16 0.16"
  cor_val = c(cor_val,get_cor_two_groups(sce,group,grp2))
  
  group = "0.33 0.33 0.33"
  grp2=c("0.68 0.16 0.16","0.16 0.68 0.16","0.16 0.16 0.68")
  cor_val = c(cor_val,get_cor_two_groups(sce,group,grp2))

  return(mean(cor_val))
}

hivar_method = list(scran_hi = scran_high_var)

norm_evaluation <- list(
  silhouette_mean=silhouette_pca,
  silhouette_noPCA=silhouette_raw,
  unwanted_variation_corr=fac_corr_u,
  wanted_variation_corr=fac_corr_w
  
)

norm_evaluation_RNAmix = list(
  nearest_group=cor_RNAmix_between,
  within_group=cor_RNAmix_within
)



####
get_method_time = function(sce){
  typ="imputation"
  if(!("logcounts" %in% assayNames(sce))){
    return(NA)
  }
  if (!("running_time" %in% names(metadata(sce)))){
    return(NA)
  }else if(nrow( metadata(sce)$running_time[metadata(sce)$running_time$method_type==typ,])==0){
    return(NA)
  }
  tmp = metadata(sce)$running_time
  tmp = tmp[tmp$method_type==typ,"time"]
  return(tmp)
}

get_total_time = function(sce){
  typ="imputation"
  if(!("logcounts" %in% assayNames(sce))){
    return(NA)
  }
  if (!("running_time" %in% names(metadata(sce)))){
    return(NA)
  }else if(nrow( metadata(sce)$running_time[metadata(sce)$running_time$method_type==typ,])==0){
    return(NA)
  }
  tmp = metadata(sce)$running_time
  return(sum(tmp$time))
}

get_cell_number = function(sce){
  return(ncol(sce))
}

regress_running_time = function(res){
  method_v = c()
  scale_coeff = c()
  for (i in unique(res$impute_method)){
    tmp = res[res$impute_method==i,]
    tmp$method_time = log2(tmp$method_time+1)
    tmp$cell_number = log2(tmp$cell_number+1) 
    fit = lm(method_time~cell_number,data=tmp)
    print(tmp)
    sm = summary(fit)
    print(sm)
    scale_coeff = c(scale_coeff,sm$coefficients[2,1])
    method_v = c(method_v, i)
  }
  return(data.frame(method=method_v,time_coeff=scale_coeff,stringsAsFactors = FALSE))
}

timing_method = list(method_time=get_method_time,
                     total_time=get_total_time,
                     cell_number=get_cell_number)

####


after_imputation <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/sc_all_after_imputation.Rds")
after_imputation = after_imputation[unlist(lapply(after_imputation$result, function(x){"logcounts" %in% assayNames(x)})),]
after_imputation = after_imputation %>% apply_methods(hivar_method)
after_imputation = after_imputation[,!(colnames(after_imputation) %in% c("hivar_method"))]
res = after_imputation %>% apply_methods(norm_evaluation)
res$cell_conditions="single_cell"
clu_time = after_imputation %>% apply_methods(timing_method)

after_imputation <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/cellmix_all_after_imputation.Rds")
after_imputation = after_imputation[unlist(lapply(after_imputation$result, function(x){"logcounts" %in% assayNames(x)})),]
after_imputation = after_imputation %>% apply_methods(hivar_method)
after_imputation = after_imputation[,!(colnames(after_imputation) %in% c("hivar_method"))]
tmp = after_imputation %>% apply_methods(norm_evaluation)
tmp$cell_conditions="cellmix"
res = rbind(res,tmp)
tmp = after_imputation %>% apply_methods(timing_method)
clu_time = rbind(clu_time,tmp)

after_imputation <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/RNAmix_all_after_imputation.Rds")
after_imputation = after_imputation[unlist(lapply(after_imputation$result, function(x){"logcounts" %in% assayNames(x)})),]
after_imputation = after_imputation %>% apply_methods(hivar_method)
after_imputation = after_imputation[,!(colnames(after_imputation) %in% c("hivar_method"))]
RNAmix_res = after_imputation %>% apply_methods(norm_evaluation_RNAmix)
tmp = after_imputation %>% apply_methods(norm_evaluation)
tmp$cell_conditions="RNAmix"
res = rbind(res,tmp)
tmp = after_imputation %>% apply_methods(timing_method)
clu_time = rbind(clu_time,tmp)



saveRDS(clu_time,file="norm_timing_result.Rds")
saveRDS(res,file="norm_evaluation_result.Rds")
# res <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/norm_evaluation_result.Rds")

clu_time = clu_time %>% spread(timing_method,result)
clu_time = clu_time[!is.na(clu_time$method_time),]
ggplot(data=clu_time,aes(x=log2(cell_number+1),y=log2(method_time+1),col=impute_method))+geom_point()+theme_bw()
time_coeff = regress_running_time(clu_time)




ggplot(data=res[res$norm_evaluation=="silhouette_noPCA",],aes(x=norm_method,y=result,col=cell_conditions))+geom_boxplot()+theme_bw()

norm_result = function(res){
  res$norm_result=NA
  for (tp in unique(res$cell_conditions)){
    res[res$cell_conditions==tp,]$norm_result = scale(res[res$cell_conditions==tp,]$result,center=mean(res[res$cell_conditions==tp & res$norm_method=="none" & res$impute_method=="no_imputation",]$result))
  }
  return(res)
}

res$impute_method = revalue(res$impute_method, c("no_impute"="no_imputation"))

res = res[!is.na(res$result),]
res_sil = res[res$norm_evaluation=="silhouette_noPCA",]
res_sil_norm = norm_result(res_sil)
res_corr = res[!(res$norm_evaluation=="silhouette_noPCA"),]
res_corr = res_corr %>% spread(norm_evaluation,result)
res_corr$corr_diff = res_corr$wanted_variation_corr-res_corr$unwanted_variation_corr


res_norm_sub = res_sil[res_sil$impute_method=="no_imputation",]
lv = levels(reorder(res_norm_sub$norm_method,res_norm_sub$result))
res_norm_sub$norm_method <- factor(res_norm_sub$norm_method, lv)


#pdf("norm_boxplot_noPCA_nonorm.pdf",width = 10,height = 4)
p1 = ggplot(data=res_norm_sub,aes(x=norm_method,y=result,col=cell_conditions))+
  geom_boxplot(outlier.colour = NA, position = position_dodge(width = 0.8),show.legend = TRUE)+
  geom_point(alpha=0.5,position=position_jitterdodge(dodge.width=0.8),show.legend = FALSE)+
  scale_color_brewer(palette="Set2")+
  labs(y="Silhouette width",title="Silhouette on gene expression matrix")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank())
#dev.off()


#res$impute_method = revalue(res$impute_method, c("no_impute"="no_imputation"))

res = res[!is.na(res$result),]
res_sil = res[res$norm_evaluation=="silhouette_mean",]
res_sil_norm = norm_result(res_sil)
res_corr = res[!(res$norm_evaluation=="silhouette_mean"),]
res_corr = res_corr %>% spread(norm_evaluation,result)
res_corr$corr_diff = res_corr$wanted_variation_corr-res_corr$unwanted_variation_corr


res_norm_sub = res_sil_norm[res_sil_norm$impute_method=="no_imputation",]
lv = levels(reorder(res_norm_sub$norm_method,res_norm_sub$norm_result))
res_norm_sub$norm_method <- factor(res_norm_sub$norm_method, lv)

#pdf("norm_boxplot_alldata.pdf",width = 10,height = 4)
ggplot(data=res_norm_sub,aes(x=norm_method,y=norm_result,col=norm_method))+
  geom_boxplot(outlier.colour = NA, position = position_dodge(width = 0.8),show.legend = FALSE)+
  geom_point(alpha=0.5,position=position_jitterdodge(dodge.width=0.8),show.legend = FALSE)+
  labs(y="normalized silhouette width")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank())
#dev.off()

res_norm_sub = res_sil[res_sil$impute_method=="no_imputation",]
lv = levels(reorder(res_norm_sub$norm_method,res_norm_sub$result))
res_norm_sub$norm_method <- factor(res_norm_sub$norm_method, lv)

#pdf("norm_boxplot_alldata_nonorm.pdf",width = 10,height = 4)
p2 = ggplot(data=res_norm_sub,aes(x=norm_method,y=result,col=cell_conditions))+
  geom_boxplot(outlier.colour = NA, position = position_dodge(width = 0.8),show.legend = FALSE)+
  geom_point(alpha=0.5,position=position_jitterdodge(dodge.width=0.8),show.legend = FALSE)+
  scale_color_brewer(palette="Set2")+
  labs(y="Silhouette width",title="Silhouette on PCA")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank())
#dev.off()

pdf("norm_boxplot_alldata_nonorm.pdf",width = 10,height = 8)
ggarrange(p1,p2,ncol=1,nrow=2,
          common.legend=TRUE)
dev.off()

ggplot(data=res_sil_norm[,],aes(x=impute_method,y=norm_result,col=norm_method))+geom_boxplot()+theme_bw()


ggplot(data=res_sil[res_sil$cell_conditions=="RNAmix",],aes(x=impute_method,y=result,col=norm_method))+geom_boxplot()+theme_bw()

ggplot(data=res_sil[res_sil$cell_conditions=="cellmix",],aes(x=impute_method,y=result,col=norm_method))+geom_boxplot()+theme_bw()

ggplot(data=res_sil[res_sil$cell_conditions=="single_cell",],aes(x=impute_method,y=result,col=norm_method))+geom_boxplot()+theme_bw()

####
ggplot(data=res_corr,aes(x=norm_method,y=corr_diff,col=impute_method))+geom_boxplot()+theme_bw()
####

RNAmix_res$norm_method <- factor(RNAmix_res$norm_method, lv)
RNAmix_res$impute_method = revalue(RNAmix_res$impute_method, c("no_impute"="no_imputation"))
pdf("within_group_corr_RNAmix_imputation.pdf",width = 13,height = 5)
ggplot(data=RNAmix_res[RNAmix_res$norm_evaluation_RNAmix=="within_group",],aes(x=impute_method,y=result,col=norm_method))+
  geom_boxplot(outlier.colour = NA, position = position_dodge(width = 0.8))+
  geom_point(alpha=0.5,position=position_jitterdodge(dodge.width=0.8))+
  labs(y="within-group correlation",col="normalization method")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank())
dev.off()

ix=8
plotPCA(after_imputation$result[[ix]],colour_by="group",ncomponents = 3)

#sce = after_imputation$result[[5]]
#library(SAVER)
#assay(sce,"saver_log") = saver(logcounts(after_imputation$result[[5]]), ncores=1, size.factor=1, estimates.only = TRUE)
#assay(sce,"saver_no_log")  = log2(saver(2^(logcounts(after_imputation$result[[5]]))-1, ncores=1, size.factor=1, estimates.only = TRUE)+1)
#
#saveRDS(sce,file="saver_test_sce.Rds")

fit = lm(result~cell_conditions+norm_method+impute_method,data=res_norm)


library(pheatmap)
library(RColorBrewer)
mRNA_col = brewer.pal(4, "RdPu")

sce = after_imputation[after_imputation$norm_method=="none" & after_imputation$impute_method=="knn_smooth2",]$result[[1]]
sce$mRNA_amount = as.factor(as.character( sce$mRNA_amount))
sce = sce[,sce$group %in% c("0 0 1","1 0 0")]
sce$group[sce$group=="0 0 1"] = "HCC827: 100%"
sce$group[sce$group=="1 0 0"] = "H2228: 100%"
sce = scran_high_var(sce,topn=500)
lc = logcounts(sce)[rowData(sce)$hi_var,]
ccor = cor(lc)
pdf("hm_knn_smooth2_with_legend.pdf",width=7,height = 5)
pheatmap(ccor,annotation_col=as.data.frame(colData(sce)[,c("group","mRNA_amount")]),annotation_colors=list(mRNA_amount=c("3.75"=mRNA_col[1],
                                                                                                                               "7.5"=mRNA_col[2],
                                                                                                                               "15"=mRNA_col[3],
                                                                                                                               "30"=mRNA_col[4]),
                                                                                                           group=c("H2228: 100%"="#6A8ECF",
                                                                                                             "HCC827: 100%"="#87BA66")),
         show_rownames=FALSE,
         show_colnames=FALSE,
         border_color = NA,
         treeheight_row = 0, treeheight_col = 0)
dev.off()

pdf("hm_knn_smooth2.pdf",width=5,height = 5)
pheatmap(ccor,annotation_col=as.data.frame(colData(sce)[,c("group","mRNA_amount")]),annotation_colors=list(mRNA_amount=c("3.75"=mRNA_col[1],
                                                                                                                         "7.5"=mRNA_col[2],
                                                                                                                         "15"=mRNA_col[3],
                                                                                                                         "30"=mRNA_col[4]),
                                                                                                           group=c("H2228: 100%"="#6A8ECF",
                                                                                                                   "HCC827: 100%"="#87BA66")),
         show_rownames=FALSE,
         show_colnames=FALSE,
         border_color = NA,
         treeheight_row = 0, treeheight_col = 0,annotation_legend=FALSE)
dev.off()



sce = after_imputation[after_imputation$norm_method=="Linnorm" & after_imputation$impute_method=="SAVER",]$result[[1]]
sce = sce[,!(colnames(sce)=="L19")]
sce = runPCA(sce)
col <- rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
dr = reducedDim(sce,"PCA")
pdf("scatter_Linnorm_SAVER.pdf",width=5,height = 5)
ggplot(data=as.data.frame(reducedDim(sce,"PCA")),aes(x=PC1,y=PC2,col=sce$group,size=sce$mRNA_amount))+
  geom_point(show.legend = F)+
  scale_size_continuous(range = c(1,4))+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=3),"%)"))+
  scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()


sce$mRNA_amount = as.factor(as.character( sce$mRNA_amount))
sce = sce[,sce$group %in% c("0 0 1","1 0 0")]
sce$group[sce$group=="0 0 1"] = "HCC827: 100%"
sce$group[sce$group=="1 0 0"] = "H2228: 100%"
sce = scran_high_var(sce,topn=500)
lc = logcounts(sce)[rowData(sce)$hi_var,]
ccor = cor(lc)
pdf("hm_Linnorm_SAVER.pdf",width=5,height = 5)
pheatmap(ccor,annotation_col=as.data.frame(colData(sce)[,c("group","mRNA_amount")]),annotation_colors=list(mRNA_amount=c("3.75"=mRNA_col[1],
                                                                                                                         "7.5"=mRNA_col[2],
                                                                                                                         "15"=mRNA_col[3],
                                                                                                                         "30"=mRNA_col[4]),
                                                                                                           group=c("H2228: 100%"="#6A8ECF",
                                                                                                                   "HCC827: 100%"="#87BA66")),
         show_rownames=FALSE,
         show_colnames=FALSE,
         border_color = NA,
         treeheight_row = 0, treeheight_col = 0,annotation_legend=FALSE)
dev.off()





sce = after_imputation[after_imputation$norm_method=="Linnorm" & after_imputation$impute_method=="DrImpute",]$result[[1]]
sce$mRNA_amount = as.factor(as.character( sce$mRNA_amount))
sce = sce[,sce$group %in% c("0 0 1","1 0 0")]
sce$group[sce$group=="0 0 1"] = "HCC827: 100%"
sce$group[sce$group=="1 0 0"] = "H2228: 100%"
sce = scran_high_var(sce,topn=500)
lc = logcounts(sce)[rowData(sce)$hi_var,]
ccor = cor(lc)
pdf("hm_Linnorm_DrImpute.pdf",width=5,height = 5)
pheatmap(ccor,annotation_col=as.data.frame(colData(sce)[,c("group","mRNA_amount")]),annotation_colors=list(mRNA_amount=c("3.75"=mRNA_col[1],
                                                                                                                         "7.5"=mRNA_col[2],
                                                                                                                         "15"=mRNA_col[3],
                                                                                                                         "30"=mRNA_col[4]),
                                                                                                           group=c("H2228: 100%"="#6A8ECF",
                                                                                                                   "HCC827: 100%"="#87BA66")),
         show_rownames=FALSE,
         show_colnames=FALSE,
         border_color = NA,
         treeheight_row = 0, treeheight_col = 0,annotation_legend=FALSE)
dev.off()


sce = after_imputation[after_imputation$norm_method=="Linnorm" & after_imputation$impute_method=="no_impute",]$result[[1]]

sce = sce[,!(colnames(sce)=="L19")]
sce = runPCA(sce)
col <- rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
dr = reducedDim(sce,"PCA")
pdf("scatter_Linnorm_no_imputation.pdf",width=5,height = 5)
ggplot(data=as.data.frame(reducedDim(sce,"PCA")),aes(x=PC1,y=PC2,col=sce$group,size=sce$mRNA_amount))+
  geom_point(show.legend = F)+
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=3),"%)"))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

sce$mRNA_amount = as.factor(as.character( sce$mRNA_amount))
sce = sce[,sce$group %in% c("0 0 1","1 0 0")]
sce$group[sce$group=="0 0 1"] = "HCC827: 100%"
sce$group[sce$group=="1 0 0"] = "H2228: 100%"
sce = scran_high_var(sce,topn=500)
lc = logcounts(sce)[rowData(sce)$hi_var,]
ccor = cor(lc)
pdf("hm_Linnorm_no_imputation.pdf",width=5,height = 5)
pheatmap(ccor,annotation_col=as.data.frame(colData(sce)[,c("group","mRNA_amount")]),annotation_colors=list(mRNA_amount=c("3.75"=mRNA_col[1],
                                                                                                                         "7.5"=mRNA_col[2],
                                                                                                                         "15"=mRNA_col[3],
                                                                                                                         "30"=mRNA_col[4]),
                                                                                                           group=c("H2228: 100%"="#6A8ECF",
                                                                                                                   "HCC827: 100%"="#87BA66")),
         show_rownames=FALSE,
         show_colnames=FALSE,
         border_color = NA,
         treeheight_row = 0, treeheight_col = 0,annotation_legend=FALSE)
dev.off()


sce = after_imputation[after_imputation$norm_method=="scran" & after_imputation$impute_method=="knn_smooth2",]$result[[1]]
sce = sce[,!(colnames(sce)=="L19")]
sce = runPCA(sce)
col <- rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
dr = reducedDim(sce,"PCA")
pdf("scatter_scran_knn_smooth2.pdf",width=5,height = 5)
ggplot(data=as.data.frame(reducedDim(sce,"PCA")),aes(x=PC1,y=PC2,col=sce$group,size=sce$mRNA_amount))+
  geom_point(show.legend = F)+
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=3),"%)"))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()



sce = after_imputation[after_imputation$norm_method=="scone" & after_imputation$impute_method=="SAVER",]$result[[1]]
sce = sce[,!(colnames(sce)=="L19")]
sce = runPCA(sce)
col <- rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
dr = reducedDim(sce,"PCA")

pdf("scatter_scone_SAVER.pdf",width=5,height = 5)
ggplot(data=as.data.frame(reducedDim(sce,"PCA")),aes(x=PC1,y=PC2,col=sce$group,size=sce$mRNA_amount))+
  geom_point(show.legend = F)+
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=3),"%)"))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()


sce = after_imputation[after_imputation$norm_method=="none" & after_imputation$impute_method=="knn_smooth2",]$result[[1]]
sce = sce[,!(colnames(sce)=="L19")]
sce = runPCA(sce)
col <- rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
dr = reducedDim(sce,"PCA")

pdf("scatter_none_knn_smooth2.pdf",width=5,height = 5)
ggplot(data=as.data.frame(reducedDim(sce,"PCA")),aes(x=PC1,y=PC2,col=sce$group,size=sce$mRNA_amount))+
  geom_point(show.legend = F)+
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=3),"%)"))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()



sce = after_imputation[after_imputation$norm_method=="scran" & after_imputation$impute_method=="DrImpute",]$result[[1]]
sce = sce[,!(colnames(sce)=="L19")]
sce = runPCA(sce)
col <- rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
dr = reducedDim(sce,"PCA")

pdf("scatter_scran_DrImpute.pdf",width=5,height = 5)
ggplot(data=as.data.frame(reducedDim(sce,"PCA")),aes(x=PC1,y=PC2,col=sce$group,size=sce$mRNA_amount))+
  geom_point(show.legend = F)+
  scale_size_continuous(range = c(1,4))+
  scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=3),"%)"))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()



