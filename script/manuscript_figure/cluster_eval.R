setwd("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit")

library(CellBench)
library(ggplot2)
library(scater)
library(scran)
library(tidyr)
library(plyr)
library(cluster)
library(mclust)
library(ggpubr)
library(RColorBrewer)
library(umap)

ARI_matric = function(sce){
  if(!("clustering_res" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    ari_val = adjustedRandIndex(sce$group, sce$clustering_res)
  }else{
    ari_val = adjustedRandIndex(sce$cell_line_demuxlet, sce$clustering_res)
  }
  
  return(ari_val)
}

cluster_number = function(sce){
  if(!("clustering_res" %in% colnames(colData(sce)))){
    return(NA)
  }
  return(length(table(sce$clustering_res)))
}

cal_entropy=function(x){
  freqs <- table(x)/length(x)
  freqs = freqs[freqs>0]
  return(-sum(freqs * log(freqs)))
}

get_cluster_purity=function(sce){
  if(!("clustering_res" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    return(mean(unlist(lapply(unique(sce$group),function(x){cal_entropy(sce$clustering_res[sce$group==x])}))))
  }else{
    return(mean(unlist(lapply(unique(sce$cell_line_demuxlet),function(x){cal_entropy(sce$clustering_res[sce$cell_line_demuxlet==x])}))))
  }
  
}

get_cluster_accuracy=function(sce){
  if(!("clustering_res" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    return(mean(unlist(lapply(unique(sce$clustering_res),function(x){cal_entropy(sce$group[sce$clustering_res==x])}))))
  }else{
    return(mean(unlist(lapply(unique(sce$clustering_res),function(x){cal_entropy(sce$cell_line_demuxlet[sce$clustering_res==x])}))))
  }
  
}


clustering_evaluation <- list(
  ARI=ARI_matric,
  cluster_purity=get_cluster_purity,
  cluster_accuracy=get_cluster_accuracy,
  cluster_number=cluster_number
)



get_method_time = function(sce){
  typ="clustering"
  if(!("clustering_res" %in% colnames(colData(sce)))){
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
  typ="clustering"
  if(!("clustering_res" %in% colnames(colData(sce)))){
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
  for (i in unique(res$clustering_method)){
    tmp = res[res$clustering_method==i,]
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

# clu_time = after_clustering %>% apply_methods(timing_method)
# clu_time = clu_time %>% spread(timing_method,result)
# clu_time = clu_time[!is.na(clu_time$method_time),]
# ggplot(data=clu_time,aes(x=log2(cell_number+1),y=log2(method_time+1),col=clustering_method))+geom_point()+theme_bw()
# time_coeff = regress_running_time(clu_time)






after_clustering <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/cellmix_all_after_clustering.Rds")
res = after_clustering %>% apply_methods(clustering_evaluation)
clu_time = after_clustering %>% apply_methods(timing_method)
res$cell_conditions="cellmix"
after_clustering <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/cellmix_after_clustering_mac.Rds")
tmp = after_clustering %>% apply_methods(clustering_evaluation)
tmp$cell_conditions="cellmix"
res = rbind(res,tmp)
tmp = after_clustering %>% apply_methods(timing_method)
clu_time = rbind(clu_time,tmp)


after_clustering <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/sc_all_after_clustering.Rds")
tmp = after_clustering %>% apply_methods(clustering_evaluation)
tmp$cell_conditions="single_cell_3cl"
tmp$cell_conditions[grepl("5cl",tmp$data)]="single_cell_5cl"
res = rbind(res,tmp)
tmp = after_clustering %>% apply_methods(timing_method)
clu_time = rbind(clu_time,tmp)



after_clustering <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/RNAmix_all_after_clustering.Rds")
tmp = after_clustering %>% apply_methods(clustering_evaluation)
tmp$cell_conditions="RNAmix"
res = rbind(res,tmp)
tmp = after_clustering %>% apply_methods(timing_method)
clu_time = rbind(clu_time,tmp)

after_clustering <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/RNAmix_after_clustering_mac.Rds")
tmp = after_clustering %>% apply_methods(clustering_evaluation)
tmp$cell_conditions="RNAmix"
res = rbind(res,tmp)
tmp = after_clustering %>% apply_methods(timing_method)
clu_time = rbind(clu_time,tmp)


for(i in 1:30){
  after_clustering = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/from_mac/sc_after_clu_",i,".Rds"))
  tmp = after_clustering %>% apply_methods(clustering_evaluation)
  tmp$cell_conditions="single_cell_3cl"
  tmp$cell_conditions[grepl("5cl",tmp$data)]="single_cell_5cl"
  res = rbind(res,tmp)
  tmp = after_clustering %>% apply_methods(timing_method)
  clu_time = rbind(clu_time,tmp)
}


saveRDS(res,file="cluster_evaluation_result.Rds")
saveRDS(clu_time,file="cluster_timing_result.Rds")



clu_time = clu_time %>% spread(timing_method,result)
clu_time = clu_time[!is.na(clu_time$method_time),]
ggplot(data=clu_time,aes(x=log2(cell_number+1),y=log2(method_time+1),col=clustering_method))+geom_point()+theme_bw()
time_coeff = regress_running_time(clu_time)

res_spread = res %>% spread(clustering_evaluation,result)
res_spread = res_spread[!is.na(res_spread$cluster_purity),]
plot(res_spread$cluster_purity,res_spread$cluster_number)


DE_methods = function(res_spread){
  res_method = c()
  res_m_type = c()
  res_coeff = c()
  res_pval = c()
  for(i in unique(as.character(res_spread$impute_method))){
    res_spread$the_method = "NO"
    res_spread$the_method[res_spread$impute_method==i] = "YES"
    fit = lm(ARI ~ the_method+cell_conditions,data=res_spread)
    sm = summary(fit)
    res_method = c(res_method,i)
    res_m_type = c(res_m_type,"impute_method")
    res_coeff = c(res_coeff,sm$coefficients[2,1])
    res_pval= c(res_pval,sm$coefficients[2,4])
  }
  for(i in unique(as.character(res_spread$norm_method))){
    res_spread$the_method = "NO"
    res_spread$the_method[res_spread$norm_method==i] = "YES"
    fit = lm(ARI ~ the_method+cell_conditions,data=res_spread)
    sm = summary(fit)
    res_method = c(res_method,i)
    res_m_type = c(res_m_type,"norm_method")
    res_coeff = c(res_coeff,sm$coefficients[2,1])
    res_pval= c(res_pval,sm$coefficients[2,4])
  }
  for(i in unique(as.character(res_spread$clustering_method))){
    res_spread$the_method = "NO"
    res_spread$the_method[res_spread$clustering_method==i] = "YES"
    fit = lm(ARI ~ the_method+cell_conditions,data=res_spread)
    sm = summary(fit)
    res_method = c(res_method,i)
    res_m_type = c(res_m_type,"clustering_method")
    res_coeff = c(res_coeff,sm$coefficients[2,1])
    res_pval= c(res_pval,sm$coefficients[2,4])
  }
  res_df = data.frame(method=res_method,
                      method_type=res_m_type,
                      coefficient=res_coeff,
                      p_value=res_pval,
                      stringsAsFactors = FALSE)
  return(res_df)
}

DE_methods_clu_num = function(res_spread){
  res_method = c()
  res_m_type = c()
  res_coeff = c()
  res_pval = c()
  for(i in unique(as.character(res_spread$impute_method))){
    res_spread$the_method = "NO"
    res_spread$the_method[res_spread$impute_method==i] = "YES"
    fit = lm(cluster_number ~ the_method+cell_conditions,data=res_spread)
    sm = summary(fit)
    res_method = c(res_method,i)
    res_m_type = c(res_m_type,"impute_method")
    res_coeff = c(res_coeff,sm$coefficients[2,1])
    res_pval= c(res_pval,sm$coefficients[2,4])
  }
  for(i in unique(as.character(res_spread$norm_method))){
    res_spread$the_method = "NO"
    res_spread$the_method[res_spread$norm_method==i] = "YES"
    fit = lm(cluster_number ~ the_method+cell_conditions,data=res_spread)
    sm = summary(fit)
    res_method = c(res_method,i)
    res_m_type = c(res_m_type,"norm_method")
    res_coeff = c(res_coeff,sm$coefficients[2,1])
    res_pval= c(res_pval,sm$coefficients[2,4])
  }
  for(i in unique(as.character(res_spread$clustering_method))){
    res_spread$the_method = "NO"
    res_spread$the_method[res_spread$clustering_method==i] = "YES"
    fit = lm(cluster_number ~ the_method+cell_conditions,data=res_spread)
    sm = summary(fit)
    res_method = c(res_method,i)
    res_m_type = c(res_m_type,"clustering_method")
    res_coeff = c(res_coeff,sm$coefficients[2,1])
    res_pval= c(res_pval,sm$coefficients[2,4])
  }
  res_df = data.frame(method=res_method,
                      method_type=res_m_type,
                      coefficient=res_coeff,
                      p_value=res_pval,
                      stringsAsFactors = FALSE)
  return(res_df)
}

ggplot(data=res_spread,aes(x=clustering_method,y=cluster_number,col=cell_conditions))+
  geom_boxplot()+
  theme_bw()
head(res_spread[order(res_spread$cluster_number,decreasing = TRUE),])

fit = lm(cluster_number~clustering_method+norm_method+impute_method+cell_conditions,data=res_spread)
summary(fit)

fit = lm(ARI~clustering_method+norm_method+impute_method+cell_conditions,data=res_spread)
summary(fit)


ggplot(data=res_spread,aes(x=clustering_method,y=ARI,col=cell_conditions))+
  geom_boxplot()+
  theme_bw()

DE_res = DE_methods(res_spread)
DE_res = DE_res[order(DE_res$coefficient),]
DE_res = rbind(DE_res[DE_res$method_type=="norm_method",],DE_res[DE_res$method_type=="impute_method",],DE_res[DE_res$method_type=="clustering_method",])
DE_res$method_type[DE_res$method_type=="clustering_method"] = "clustering"
DE_res$method_type[DE_res$method_type=="norm_method"] = "normalization"
DE_res$method_type[DE_res$method_type=="impute_method"] = "imputation"
DE_res$method = factor(DE_res$method, levels = unique(DE_res$method))
pdf("coeff_ARI_clu.pdf",width = 11,height = 5)
ggplot(data=DE_res,aes(x=method,y=coefficient,fill=method_type))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette="Set2")+
  ggtitle("ARI")+
  labs(fill="")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()


DE_res = DE_methods_clu_num(res_spread)
DE_res = DE_res[order(DE_res$coefficient),]
DE_res = rbind(DE_res[DE_res$method_type=="norm_method",],DE_res[DE_res$method_type=="impute_method",],DE_res[DE_res$method_type=="clustering_method",])
DE_res$method_type[DE_res$method_type=="clustering_method"] = "clustering"
DE_res$method_type[DE_res$method_type=="norm_method"] = "normalization"
DE_res$method_type[DE_res$method_type=="impute_method"] = "imputation"
DE_res$method = factor(DE_res$method, levels = unique(DE_res$method))
pdf("coeff_cluster_number_clu.pdf",width = 11,height = 5)
ggplot(data=DE_res,aes(x=method,y=coefficient,fill=method_type))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette="Set2")+
  ggtitle("number of clusters")+
  labs(fill="")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()


res_ent = res[(res$clustering_evaluation %in% c("cluster_purity","cluster_accuracy")),]
res_ent = res_ent[!is.na(res_ent$result),]
res_ari = res[res$clustering_evaluation %in% c("ARI", "cluster_number"),]
res_ari = res_ari[!is.na(res_ari$result),]

ggplot(data=res_ent[res_ent$cell_conditions=="RNAmix",],aes(x=clustering_method,y=result,col=clustering_evaluation))+
  geom_boxplot()+
  theme_bw()

ARi_sel_top = function(res,topn=3){
  top_res = NULL
  for (i in unique(res$clustering_method)){
    tmp = res[res$clustering_method==i,]
    tmp = tmp[order(tmp$ARI,decreasing = TRUE)[1:min(topn,nrow(tmp))],]
    if(is.null(top_res)){
      top_res=tmp
    }else{
      top_res = rbind(tmp,top_res)
    }
  }
  return(top_res)
}


tmp = res_ari[res_ari$cell_conditions=="cellmix",]
tmp = tmp %>% spread(clustering_evaluation,result)
tmp = ARi_sel_top(tmp)
tmp$norm_imputation = paste(tmp$norm_method,tmp$impute_method,sep="-")
ggplot(data=tmp,aes(x=cluster_number,y=ARI,shape=clustering_method,col=clustering_method))+
  geom_point(size=5)+
  labs(y="ARI")+
  scale_shape_manual(values=1:nlevels(tmp$clustering_method))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank())

tmp = res_ari[res_ari$cell_conditions=="RNAmix",]
tmp = tmp %>% spread(clustering_evaluation,result)
tmp = ARi_sel_top(tmp)
tmp$norm_imputation = paste(tmp$norm_method,tmp$impute_method,sep="-")
ggplot(data=tmp,aes(x=cluster_number,y=ARI,shape=clustering_method,col=clustering_method))+
  geom_point(size=5)+
  labs(y="ARI")+
  scale_shape_manual(values=1:nlevels(tmp$clustering_method))+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank())

############### 

entropy_sel_top = function(res,topn=2){
  top_res = NULL
  for (i in unique(res$clustering_method)){
    tmp = res[res$clustering_method==i,]
    tmp = tmp[order(tmp$entropy_sum,decreasing = FALSE)[1:min(topn,nrow(tmp))],]
    if(is.null(top_res)){
      top_res=tmp
    }else{
      top_res = rbind(tmp,top_res)
    }
  }
  return(top_res)
}

entropy_top_per_data = function(res){
  top_res = NULL
  for (i in unique(res$clustering_method)){
    for(j in unique(res$data)){
      tmp = res[res$clustering_method==i & res$data==j,]
      tmp = tmp[order(tmp$entropy_sum,decreasing = FALSE)[1],]
      if(is.null(top_res)){
        top_res=tmp
      }else{
        top_res = rbind(tmp,top_res)
      }
    }
    }

  return(top_res)
}


####### cellmix
res_ent_df = res_ent %>% spread(clustering_evaluation,result)
res_ent_df = res_ent_df[res_ent_df$cell_conditions=="cellmix",]
#res_ent_df$entropy_sum = rowSums(res_ent_df[,c("cluster_purity","cluster_accuracy")])
res_ent_df$entropy_sum = sqrt(res_ent_df$cluster_purity^2+res_ent_df$cluster_accuracy^2)
res_ent_df_top=entropy_sel_top(res_ent_df)
res_ent_df_top$norm_impute = paste(res_ent_df_top$norm_method,res_ent_df_top$impute_method,sep="_")
p1=ggplot(data=res_ent_df_top,aes(x=cluster_purity,y=cluster_accuracy,col=clustering_method ,shape=clustering_method ,label=norm_impute))+
  geom_point(size=2)+
  geom_text_repel()+
  scale_colour_brewer(palette="Dark2")+
  scale_shape_manual(values=1:nlevels(res_ent_df_top$clustering_method))+
  labs(x="entropy of cluster purity",y="entropy of cluster accuracy",col="",shape="")+
  theme_bw()+
  theme(text = element_text(size=15),legend.position="top")
p1

res_ent_df = res_ent_df %>% gather(clustering_evaluation, result, cluster_purity, cluster_accuracy)
p11 = ggplot(data=res_ent_df,aes(x=clustering_method,y=result,col=clustering_evaluation))+
  geom_boxplot(outlier.colour = NA, position = position_dodge(width = 0.8))+
  geom_point(alpha=0.5,position=position_jitterdodge(dodge.width=0.8,jitter.width=0.1))+
  labs(y="entropy",col="type")+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank())
p11

####### RNAmix
res_ent_df = res_ent %>% spread(clustering_evaluation,result)
res_ent_df = res_ent_df[res_ent_df$cell_conditions=="RNAmix",]
res_ent_df$entropy_sum = rowSums(res_ent_df[,c("cluster_purity","cluster_accuracy")])
res_ent_df_top=entropy_sel_top(res_ent_df,topn=2)
res_ent_df_top$norm_impute = paste(res_ent_df_top$norm_method,res_ent_df_top$impute_method,sep="_")
p2=ggplot(data=res_ent_df_top,aes(x=cluster_purity,y=cluster_accuracy,col=clustering_method ,shape=clustering_method ,label=norm_impute))+
  geom_point(size=2)+
  geom_text_repel()+
  scale_colour_brewer(palette="Dark2")+
  scale_shape_manual(values=1:nlevels(res_ent_df_top$clustering_method))+
  labs(x="entropy of cluster purity",y="entropy of cluster accuracy",col="",shape="")+
  theme_bw()+
  theme(text = element_text(size=15),legend.position="top")
p2

res_ent_df = res_ent_df %>% gather(clustering_evaluation, result, cluster_purity, cluster_accuracy)
p21 = ggplot(data=res_ent_df,aes(x=clustering_method,y=result,col=clustering_evaluation))+
  geom_boxplot(outlier.colour = NA, position = position_dodge(width = 0.8))+
  geom_point(alpha=0.5,position=position_jitterdodge(dodge.width=0.8,jitter.width=0.1))+
  labs(y="entropy",col="type")+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank())
p21

####### single_cell_5cl
res_ent_df = res_ent %>% spread(clustering_evaluation,result)
res_ent_df = res_ent_df[res_ent_df$cell_conditions=="single_cell_5cl",]
res_ent_df$entropy_sum = rowSums(res_ent_df[,c("cluster_purity","cluster_accuracy")])
res_ent_df_top=entropy_sel_top(res_ent_df,topn=2)
res_ent_df_top$norm_impute = paste(res_ent_df_top$norm_method,res_ent_df_top$impute_method,sep="_")
p4=ggplot(data=res_ent_df_top,aes(x=cluster_purity,y=cluster_accuracy,col=clustering_method ,shape=clustering_method ,label=norm_impute))+
  geom_point(size=2)+
  geom_text_repel()+
  scale_colour_brewer(palette="Dark2")+
  scale_shape_manual(values=1:nlevels(res_ent_df_top$clustering_method))+
  labs(x="entropy of cluster purity",y="entropy of cluster accuracy",col="",shape="")+
  theme_bw()+
  theme(text = element_text(size=15),legend.position="top")
p4

res_ent_df = res_ent_df %>% gather(clustering_evaluation, result, cluster_purity, cluster_accuracy)
p41 = ggplot(data=res_ent_df,aes(x=clustering_method,y=result,col=clustering_evaluation))+
  geom_boxplot(outlier.colour = NA, position = position_dodge(width = 0.8))+
  geom_point(alpha=0.5,position=position_jitterdodge(dodge.width=0.8,jitter.width=0.1))+
  labs(y="entropy",col="normalization method")+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank())
p41

####### single_cell_3cl
res_ent_df = res_ent %>% spread(clustering_evaluation,result)
res_ent_df = res_ent_df[res_ent_df$cell_conditions=="single_cell_3cl",]
res_ent_df$entropy_sum = rowSums(res_ent_df[,c("cluster_purity","cluster_accuracy")])
res_ent_df_top=entropy_sel_top(res_ent_df,topn=2)
res_ent_df_top$norm_impute = paste(res_ent_df_top$norm_method,res_ent_df_top$impute_method,sep="_")
p3=ggplot(data=res_ent_df_top,aes(x=cluster_purity,y=cluster_accuracy,col=clustering_method ,shape=clustering_method ,label=norm_impute))+
  geom_point(size=2)+
  geom_text_repel()+
  scale_colour_brewer(palette="Dark2")+
  scale_shape_manual(values=1:nlevels(res_ent_df_top$clustering_method))+
  labs(x="entropy of cluster purity",y="entropy of cluster accuracy",col="",shape="")+
  theme_bw()+
  theme(text = element_text(size=15),legend.position="top")
p3

res_ent_df = res_ent_df %>% gather(clustering_evaluation, result, cluster_purity, cluster_accuracy)
ggplot(data=res_ent_df,aes(x=clustering_method,y=result,col=clustering_evaluation))+
  geom_boxplot(outlier.colour = NA, position = position_dodge(width = 0.8))+
  geom_point(alpha=0.5,position=position_jitterdodge(dodge.width=0.8,jitter.width=0.1))+
  labs(y="entropy",col="normalization method")+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank())

pdf("clu_summary_all_ent.pdf",width = 10,height = 10)
ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2,common.legend  =TRUE,align="hv")
dev.off()






#scatter plots
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#### RNAmix
sce_after_imputation <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/RNAmix_all_after_imputation.Rds")
after_clustering <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/RNAmix_all_after_clustering.Rds")
sce = 
  sce_after_imputation[sce_after_imputation$norm_method=="Linnorm" & sce_after_imputation$impute_method=="SAVER" & sce_after_imputation$data=="RNAmix_Sortseq",]$result[[1]]

sce = runPCA(sce)
sce = runTSNE(sce)

res_sel = after_clustering[after_clustering$norm_method=="Linnorm" & after_clustering$impute_method=="SAVER" & after_clustering$data=="RNAmix_Sortseq" & after_clustering$clustering_method=="RaceID3",]$result[[1]]
dr = reducedDim(sce,"PCA")
pdf("PCA_Linnorm_SAVER_RaceID3_RNAmix.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=as.factor(res_sel$clustering_res)))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  scale_color_brewer(palette="Set1")+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=2),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=2),"%)"))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()
dr = reducedDim(sce,"TSNE")
pdf("tSNE_Linnorm_SAVER_RaceID3_RNAmix.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=as.factor(res_sel$clustering_res)))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  scale_color_brewer(palette="Set1")+
  labs(x="Dim1",y="Dim2")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

###cellmix

sce_after_imputation <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/cellmix_all_after_imputation.Rds")
after_clustering <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/cellmix_after_clustering_mac.Rds")
sce = 
  sce_after_imputation[sce_after_imputation$norm_method=="Linnorm" & sce_after_imputation$impute_method=="SAVER" & sce_after_imputation$data=="cellmix4",]$result[[1]]

sce = runPCA(sce)
sce = runTSNE(sce)

res_sel = after_clustering[after_clustering$norm_method=="Linnorm" & after_clustering$impute_method=="SAVER" & after_clustering$data=="cellmix4" & after_clustering$clustering_method=="clusterExperiment",]$result[[1]]
dr = reducedDim(sce,"PCA")
pdf("PCA_Linnorm_SAVER_clusterExperiment_cellmix4.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=as.factor(res_sel$clustering_res)))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  scale_color_manual(values = getPalette(length(unique(res_sel$clustering_res))))+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=2),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=2),"%)"))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()
dr = reducedDim(sce,"TSNE")
pdf("tSNE_Linnorm_SAVER_clusterExperiment_cellmix4.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=as.factor(res_sel$clustering_res)))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  scale_color_manual(values = getPalette(length(unique(res_sel$clustering_res))))+
  labs(x="Dim1",y="Dim2")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()


###sc_3cl

sce_after_imputation <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/sc_all_after_imputation.Rds")
after_clustering <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/sc_all_after_clustering.Rds")
sce = 
  sce_after_imputation[sce_after_imputation$norm_method=="scran" & sce_after_imputation$impute_method=="no_impute" & sce_after_imputation$data=="sc_CELseq2",]$result[[1]]

sce = runPCA(sce)
sce = runTSNE(sce)

res_sel = after_clustering[after_clustering$norm_method=="scran" & after_clustering$impute_method=="no_impute" & after_clustering$data=="sc_CELseq2" & after_clustering$clustering_method=="Seurat_0.6",]$result[[1]]
dr = reducedDim(sce,"PCA")
pdf("PCA_scran_Seurat_0.6_sc_CELseq2.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=as.factor(res_sel$clustering_res)))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  scale_color_brewer(palette="Set1")+
  #scale_color_manual(values = getPalette(length(unique(res_sel$clustering_res))))+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=2),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=2),"%)"))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()
dr = reducedDim(sce,"TSNE")
pdf("TSNE_scran_Seurat_0.6_sc_CELseq2.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=as.factor(res_sel$clustering_res)))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  scale_color_brewer(palette="Set1")+
  #scale_color_manual(values = getPalette(length(unique(res_sel$clustering_res))))+
  labs(x="Dim1",y="Dim2")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()


## sc_5cl

sce = 
  sce_after_imputation[sce_after_imputation$norm_method=="scran" & sce_after_imputation$impute_method=="no_impute" & sce_after_imputation$data=="sc_10x_5cl",]$result[[1]]
sce = scran_high_var(sce)
sce = runPCA(sce,ncomponents = 3)
umap_out = umap(t(logcounts(sce)[rowData(sce)$hi_var,]))
reducedDim(sce,"UMAP")= umap_out$layout

res_sel = after_clustering[after_clustering$norm_method=="none" & after_clustering$impute_method=="no_impute" & after_clustering$data=="sc_10x_5cl" & after_clustering$clustering_method=="Seurat_pipe",]$result[[1]]
dr = reducedDim(sce,"PCA")
pdf("PCA_scran_Seurat_0.6_sc_10x_5cl.pdf")
ggplot(data=NULL,aes(x=dr[,2],y=dr[,3],col=as.factor(res_sel$clustering_res)))+
  geom_point(size=4,alpha=1,show.legend = FALSE)+
  scale_color_manual(values = getPalette(length(unique(res_sel$clustering_res))))+
  #scale_color_manual(values = getPalette(length(unique(res_sel$clustering_res))))+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=2),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=2),"%)"))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()
dr = reducedDim(sce,"UMAP")
pdf("UMAP_scran_Seurat_0.6_sc_10x_5cl.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=as.factor(res_sel$clustering_res)))+
  geom_point(size=4,alpha=1,show.legend = FALSE)+
  scale_color_brewer(palette="Set1")+
  #scale_color_manual(values = getPalette(length(unique(res_sel$clustering_res))))+
  labs(x="Dim1",y="Dim2")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()



#### modeling:

res_ent_df = res_ent %>% spread(clustering_evaluation,result)

fit = lm(cluster_purity~clustering_method+norm_method+impute_method+cell_conditions,data=res_ent_df)
summary(fit)
anova(fit)


fit = lm(cluster_accuracy~clustering_method+norm_method+impute_method+cell_conditions,data=res_ent_df)
summary(fit)
anova(fit)
