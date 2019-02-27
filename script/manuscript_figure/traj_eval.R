#eval_trajectory
library(CellBench)
library(ggplot2)
library(scater)
library(scran)
library(tidyr)
library(plyr)
library(mclust)
library(ggpubr)
library(RColorBrewer)

setwd("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit")
get_trajectory_corr = function(sce){
  if("traj_eval" %in% names(metadata(sce))){
    return(mean(metadata(sce)$traj_eval$corr))
  }else{
    return(NA)
  }
}

get_trajectory_overlap = function(sce){
  if("traj_eval" %in% names(metadata(sce))){
    return(mean(metadata(sce)$traj_eval$overlap))
  }else{
    return(NA)
  }
}

traj_eval_method = list(traj_overlap=get_trajectory_overlap,
                        traj_corr=get_trajectory_corr)

res_top_per_data = function(res){
  top_res = NULL
  for (i in unique(res$traj_method)){
    for(j in unique(res$data)){
      tmp = res[res$traj_method==i & res$data==j,]
      tmp = tmp[order(tmp$traj_corr+tmp$traj_overlap,decreasing = TRUE)[1:2],]
      if(is.null(top_res)){
        top_res=tmp
      }else{
        top_res = rbind(tmp,top_res)
      }
    }
  }
  
  return(top_res)
}




####
get_method_time = function(sce){
  typ="trajectory"
  if(!("traj_eval" %in% names(metadata(sce)))){
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
  typ="trajectory"
  if(!("traj_eval" %in% names(metadata(sce)))){
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
  for (i in unique(res$traj_method)){
    tmp = res[res$traj_method==i,]
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


RNAmix_after_trajectory <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/RNAmix_after_trajectory.Rds")
cellmix_all_after_trajectory <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/cellmix_all_after_trajectory.Rds")
clu_time = RNAmix_after_trajectory %>% apply_methods(timing_method)
tmp = cellmix_all_after_trajectory %>% apply_methods(timing_method)
clu_time = rbind(clu_time,tmp)
saveRDS(clu_time,file="trajectory_timing_result.Rds")
clu_time = clu_time %>% spread(timing_method,result)
clu_time = clu_time[!is.na(clu_time$method_time),]
ggplot(data=clu_time,aes(x=log2(cell_number+1),y=log2(method_time+1),col=traj_method))+geom_point()+theme_bw()
time_coeff = regress_running_time(clu_time)



traj_eval_res1 = RNAmix_after_trajectory %>%  apply_methods(traj_eval_method)
traj_eval_res1$design = "RNAmix"
traj_eval_res2 = cellmix_all_after_trajectory %>%  apply_methods(traj_eval_method)
traj_eval_res2$design = "cellmix"
traj_eval_res = rbind(traj_eval_res1,traj_eval_res2)
saveRDS(traj_eval_res,file="trajectory_evaluation_result.Rds")
traj_eval_res = traj_eval_res[!is.na(traj_eval_res$result),]


traj_res_wide = traj_eval_res %>% spread(traj_eval_method,result)

traj_res_wide=traj_res_wide[!(traj_res_wide$traj_method=="Monocle2_SimplePPT"),]

traj_res_wide = res_top_per_data(traj_res_wide)

tmp = aggregate(traj_corr ~ traj_method, data=traj_res_wide,mean)
traj_res_wide$traj_method <- factor(traj_res_wide$traj_method, levels = tmp$traj_method[order(tmp$traj_corr)])

pdf("traj_corr.pdf",width = 8,height = 4)
ggplot(data=traj_res_wide,aes(x=traj_method,y=traj_corr,fill=design))+
  scale_fill_manual(values = c("#1B9E77","#D95F02"))+
  geom_boxplot(alpha=0.8)+
  theme_bw()+
  labs(x="",y="correlations")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),text = element_text(size=20))
dev.off()

pdf("traj_overlap.pdf",width = 8,height = 4)
ggplot(data=traj_res_wide,aes(x=traj_method,y=traj_overlap,fill=design))+
  scale_fill_manual(values = c("#1B9E77","#D95F02"))+
  geom_boxplot(alpha=0.8)+
  theme_bw()+
  labs(x="",y="overlaps")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),text = element_text(size=20))
dev.off()









####
fit = lm(result ~ 0+traj_method+impute_method+norm_method,
         data=traj_eval_res[traj_eval_res$traj_eval_method=="traj_corr",],
         na.action=na.omit)
anova(fit)

fit = lm(result ~ 0+traj_method+impute_method+norm_method,
         data=traj_eval_res[traj_eval_res$traj_eval_method=="traj_overlap",],
         na.action=na.omit)
anova(fit)

