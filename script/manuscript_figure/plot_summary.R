# plot summary

library(CellBench)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(tidyr)
setwd("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit")



#### load all result
res_int_eval_all <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/res_int_eval_all.Rds")
res_int_eval_all = res_int_eval_all[!(res_int_eval_all$data_int_method=="zinbwave"),]
cluster_evaluation_result <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/cluster_evaluation_result.Rds")
norm_evaluation_result <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/norm_evaluation_result.Rds")
trajectory_evaluation_result <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/trajectory_evaluation_result.Rds")


trajectory_timing_result <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/trajectory_timing_result.Rds")
norm_timing_result <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/norm_timing_result.Rds")
cluster_timing_result <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/cluster_timing_result.Rds")
int_timing_result <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/int_timing_result.Rds")
####

classify_speed = function(speed_vec){
  speed_discrete = rep("good",length(speed_vec))
  speed_discrete[speed_vec>log2(5*60+1)] = "fair"
  speed_discrete[speed_vec>log2(30*60+1)] = "poor"
  return(speed_discrete)
}
classify_scalability = function(sca_vec){
  sca_discrete = rep("good",length(sca_vec))
  sca_discrete[sca_vec>1] = "fair"
  sca_discrete[sca_vec>2] = "poor"
  return(sca_discrete)
}


DE_methods = function(res_spread,res_column,met_column){
  res_method = c()
  res_coeff = c()
  res_pval = c()
  for(i in unique(as.character(res_spread[,met_column]))){
    res_spread$the_method = "NO"
    res_spread$the_method[res_spread[,met_column]==i] = "YES"
    fit = lm(res_spread[,res_column] ~ res_spread$the_method+res_spread$cell_conditions)
    sm = summary(fit)
    res_method = c(res_method,i)
    res_coeff = c(res_coeff,sm$coefficients[2,1])
    res_pval= c(res_pval,sm$coefficients[2,4])
  }
  res_df = data.frame(method=res_method,
                      coefficient=res_coeff,
                      p_value=res_pval,
                      stringsAsFactors = FALSE)
  return(res_df)
}

regress_running_time = function(res,met_column){
  method_v = c()
  scale_coeff = c()
  for (i in unique(as.character(res[,met_column]))){
    print(i)
    tmp = res[res[,met_column]==i,]
    tmp$method_time = log2(tmp$method_time+1)
    tmp$cell_number = log2(tmp$cell_number+1) 
    fit = lm(method_time~cell_number,data=tmp)
    #print(tmp)
    sm = summary(fit)
    print(sm)
    scale_coeff = c(scale_coeff,sm$coefficients[2,1])
    method_v = c(method_v, i)
  }
  return(data.frame(method=method_v,time_coeff=scale_coeff,stringsAsFactors = FALSE))
}

#install.packages("xlsx")
library("xlsx")
options(java.parameters = "-Xmx2048m")

write.xlsx2(as.data.frame(norm_evaluation_result), file="Supplementary Table 4.xlsx", sheetName = "Normalization_Imputation",
            col.names = TRUE, row.names = FALSE, append = FALSE)
write.xlsx2(as.data.frame(cluster_evaluation_result), file="Supplementary Table 4.xlsx", sheetName = "Clustering",
            col.names = TRUE, row.names = FALSE, append = TRUE)
write.xlsx2(as.data.frame(trajectory_evaluation_result), file="Supplementary Table 4.xlsx", sheetName = "Trajectory",
            col.names = TRUE, row.names = FALSE, append = TRUE)
write.xlsx2(as.data.frame(res_int_eval_all), file="Supplementary Table 4.xlsx", sheetName = "Data integration",
            col.names = TRUE, row.names = FALSE, append = TRUE)

write.xlsx2(as.data.frame(norm_timing_result), file="Supplementary Table 5.xlsx", sheetName = "Normalization_Imputation",
            col.names = TRUE, row.names = FALSE, append = FALSE)
write.xlsx2(as.data.frame(cluster_timing_result), file="Supplementary Table 5.xlsx", sheetName = "Clustering",
            col.names = TRUE, row.names = FALSE, append = TRUE)
write.xlsx2(as.data.frame(trajectory_timing_result), file="Supplementary Table 5.xlsx", sheetName = "Trajectory",
            col.names = TRUE, row.names = FALSE, append = TRUE)
write.xlsx2(as.data.frame(int_timing_result), file="Supplementary Table 5.xlsx", sheetName = "Data integration",
            col.names = TRUE, row.names = FALSE, append = TRUE)


#### clustering

cluster_ARI = cluster_evaluation_result[cluster_evaluation_result$clustering_evaluation=="ARI",]
cluster_ARI = cluster_ARI[!is.na(cluster_ARI$result),]
#cluster_ARI = cluster_ARI[!(cluster_ARI$clustering_method=="Seurat_pipe"),]

cluster_sm = DE_methods(as.data.frame(cluster_ARI),res_column="result",met_column="clustering_method")

cluster_norm_sm = DE_methods(as.data.frame(cluster_ARI),res_column="result",met_column="norm_method")
cluster_impute_sm = DE_methods(as.data.frame(cluster_ARI),res_column="result",met_column="impute_method")

ARI_var= aggregate(result ~ clustering_method+cell_conditions, data=cluster_ARI, var)
ARI_var= aggregate(result ~ clustering_method, data=ARI_var, mean)

ARI_max= aggregate(result ~ clustering_method+cell_conditions, data=cluster_ARI, max)
ARI_max= aggregate(result ~ clustering_method, data=ARI_max, mean)

rownames(cluster_sm) = cluster_sm$method
rownames(ARI_var) = ARI_var$clustering_method
rownames(ARI_max) = ARI_max$clustering_method
colnames(ARI_max) = c("clustering_method","best_p")

cluster_timing_wide = cluster_timing_result %>% spread(timing_method,result)
cluster_timing_wide = cluster_timing_wide[!is.na(cluster_timing_wide$method_time),]
time_coeff = regress_running_time(as.data.frame(cluster_timing_wide),"clustering_method")
rownames(time_coeff) = time_coeff$method
time_mean = aggregate(method_time ~ clustering_method, data=cluster_timing_wide, mean)
time_mean$method_time = log2(time_mean$method_time+1)
rownames(time_mean) = time_mean$clustering_method

comb_res = cbind(cluster_sm,
                 ARI_max[rownames(cluster_sm),],
                 ARI_var[rownames(cluster_sm),],
                 time_coeff[rownames(cluster_sm),],
                 time_mean[rownames(cluster_sm),])
comb_res = comb_res[,c("best_p","coefficient","result","method_time","time_coeff")]
comb_res = comb_res[order(comb_res$best_p,decreasing = TRUE),]
colnames(comb_res) = c("Best performance", "Average performance", "Variability", "Speed", "Scalability")
comb_res_p = scale(comb_res[,1:3])
comb_res_t = comb_res[,4:5]
comb_res_t$Speed = classify_speed(comb_res_t$Speed)
comb_res_t$Scalability = classify_scalability(comb_res_t$Scalability)
anno_color = list(Speed=c("good"=RColorBrewer::brewer.pal(3,"Spectral")[1],
                          "fair"=RColorBrewer::brewer.pal(3,"Spectral")[2],
                          "poor"=RColorBrewer::brewer.pal(3,"Spectral")[3]),
                  Scalability=c("good"=RColorBrewer::brewer.pal(3,"Spectral")[1],
                                "fair"=RColorBrewer::brewer.pal(3,"Spectral")[2],
                                "poor"=RColorBrewer::brewer.pal(3,"Spectral")[3]))

# pdf("summary_clustering.pdf")
# pheatmap(comb_res_p,scale = "none",cluster_rows=FALSE,cluster_cols=FALSE,annotation_row=comb_res_t,
#          treeheight_row = 0, treeheight_col = 0,annotation_colors=anno_color,
#          legend=TRUE,show_rownames=TRUE,show_colnames=FALSE,
#          color=colorRampPalette(c(RColorBrewer::brewer.pal(3,"Spectral")[3], 
#                                   RColorBrewer::brewer.pal(3,"Spectral")[2], 
#                                   RColorBrewer::brewer.pal(3,"Spectral")[1]))(50),
#          fontsize_row=10)
# dev.off()


#### trajectory

trajectory_corr = trajectory_evaluation_result[trajectory_evaluation_result$traj_eval_method=="traj_corr",]
trajectory_corr = trajectory_corr[!is.na(trajectory_corr$result),]
trajectory_corr$cell_conditions = trajectory_corr$design

trajectory_sm = DE_methods(as.data.frame(trajectory_corr),res_column="result",met_column="traj_method")
norm_trajectory_sm = DE_methods(as.data.frame(trajectory_corr),res_column="result",met_column="norm_method")
impute_trajectory_sm = DE_methods(as.data.frame(trajectory_corr),res_column="result",met_column="impute_method")

traj_var= aggregate(result ~ traj_method+cell_conditions, data=trajectory_corr, var)
traj_var= aggregate(result ~ traj_method, data=traj_var, mean)
traj_max= aggregate(result ~ traj_method+cell_conditions, data=trajectory_corr, max)
traj_max= aggregate(result ~ traj_method, data=traj_max, mean)


rownames(trajectory_sm) = trajectory_sm$method
rownames(traj_var) = traj_var$traj_method
rownames(traj_max) = traj_max$traj_method
colnames(traj_max) = c("traj_method","best_p")

##time
trajectory_timing_wide = trajectory_timing_result %>% spread(timing_method,result)
trajectory_timing_wide = trajectory_timing_wide[!is.na(trajectory_timing_wide$method_time),]
time_coeff = regress_running_time(as.data.frame(trajectory_timing_wide),"traj_method")
rownames(time_coeff) = time_coeff$method
time_mean = aggregate(method_time ~ traj_method , data=trajectory_timing_wide, mean)
time_mean$method_time = log2(time_mean$method_time+1)
rownames(time_mean) = time_mean$traj_method 
##

comb_res = cbind(trajectory_sm,
                 traj_max[rownames(trajectory_sm),],
                 traj_var[rownames(trajectory_sm),],
                 time_coeff[rownames(trajectory_sm),],
                 time_mean[rownames(trajectory_sm),])
comb_res = comb_res[,c("best_p","coefficient","result","method_time","time_coeff")]
comb_res = comb_res[order(comb_res$best_p,decreasing = TRUE),]
colnames(comb_res) = c("Best performance", "Average performance", "Variability", "Speed", "Scalability")
comb_res_p = rbind(comb_res_p,scale(comb_res[,1:3]))
comb_res_t_tmp = comb_res[,4:5]
comb_res_t_tmp$Speed = classify_speed(comb_res_t_tmp$Speed)
comb_res_t_tmp$Scalability = classify_scalability(comb_res_t_tmp$Scalability)
comb_res_t = rbind(comb_res_t,comb_res_t_tmp)

pheatmap(comb_res_p,scale = "none",cluster_rows=FALSE,cluster_cols=FALSE,annotation_row=comb_res_t,
         treeheight_row = 0, treeheight_col = 0,annotation_colors=anno_color,
         legend=TRUE,show_rownames=TRUE,show_colnames=FALSE,
         color=colorRampPalette(c(RColorBrewer::brewer.pal(3,"Spectral")[3], 
                                  RColorBrewer::brewer.pal(3,"Spectral")[2], 
                                  RColorBrewer::brewer.pal(3,"Spectral")[1]))(50),
         gaps_row=c(7),
         fontsize_row=10)


#### data_int

data_int_sil = res_int_eval_all[res_int_eval_all$data_int_eval=="silhouette_dr",]
data_int_sil = data_int_sil[!(data_int_sil$data_int_method=="zinbwave"),]
#data_int_sil = data_int_sil[!(data_int_sil$data_int_method=="no_int"),]
data_int_sil = data_int_sil[!is.na(data_int_sil$result),]
data_int_sil$cell_conditions = data_int_sil$data

data_int_sm = DE_methods(as.data.frame(data_int_sil),res_column="result",met_column="data_int_method")
norm_data_int_sm = DE_methods(as.data.frame(data_int_sil),res_column="result",met_column="norm_method")
impute_data_int_sm = DE_methods(as.data.frame(data_int_sil),res_column="result",met_column="impute_method")

data_int_var= aggregate(result ~ data_int_method+cell_conditions, data=data_int_sil, var)
data_int_var= aggregate(result ~ data_int_method, data=data_int_var, mean)
data_int_max= aggregate(result ~ data_int_method+cell_conditions, data=data_int_sil, max)
data_int_max= aggregate(result ~ data_int_method, data=data_int_max, mean)

rownames(data_int_sm) = data_int_sm$method
rownames(data_int_var) = data_int_var$data_int_method
rownames(data_int_max) = data_int_max$data_int_method
colnames(data_int_max) = c("data_int_method","best_p")

##time
int_timing_wide = int_timing_result %>% spread(timing_method,result)
int_timing_wide = int_timing_wide[!(int_timing_wide$data_int_method=="zinbwave"),]
int_timing_wide = int_timing_wide[!is.na(int_timing_wide$method_time),]
time_coeff = regress_running_time(as.data.frame(int_timing_wide),"data_int_method")
rownames(time_coeff) = time_coeff$method
time_mean = aggregate(method_time ~ data_int_method , data=int_timing_wide, mean)
time_mean$method_time = log2(time_mean$method_time+1)
rownames(time_mean) = time_mean$data_int_method 
##

comb_res = cbind(data_int_sm,
                 data_int_max[rownames(data_int_sm),],
                 data_int_var[rownames(data_int_sm),],
                 time_coeff[rownames(data_int_sm),],
                 time_mean[rownames(data_int_sm),])
comb_res = comb_res[,c("best_p","coefficient","result","method_time","time_coeff")]
comb_res = comb_res[order(comb_res$best_p,decreasing = TRUE),]
colnames(comb_res) = c("Best performance", "Average performance", "Variability", "Speed", "Scalability")
comb_res_p = rbind(comb_res_p,scale(comb_res[,1:3]))
comb_res_t_tmp = comb_res[,4:5]
comb_res_t_tmp$Speed = classify_speed(comb_res_t_tmp$Speed)
comb_res_t_tmp$Scalability = classify_scalability(comb_res_t_tmp$Scalability)
comb_res_t = rbind(comb_res_t,comb_res_t_tmp)

comb_res_p = comb_res_p[!(rownames(comb_res_p)=="Seurat"),]

pdf("summary_clu_traj_int.pdf",width = 4,height = 8)
pheatmap(comb_res_p,scale = "none",cluster_rows=FALSE,cluster_cols=FALSE,annotation_row=comb_res_t,
         treeheight_row = 0, treeheight_col = 0,annotation_colors=anno_color,
         legend=FALSE,show_rownames=TRUE,show_colnames=TRUE,
         color=colorRampPalette(c(RColorBrewer::brewer.pal(3,"Spectral")[3], 
                                  RColorBrewer::brewer.pal(3,"Spectral")[2], 
                                  RColorBrewer::brewer.pal(3,"Spectral")[1]))(50),
         gaps_row=c(7,12),
         annotation_legend=FALSE,
         border_color="white",
         fontsize_row=10)
dev.off()
#### norm

norm_evaluation_result = norm_evaluation_result[!(norm_evaluation_result$norm_method=="none" & norm_evaluation_result$impute_method=="no_impute"),]
norm_evaluation_result = norm_evaluation_result[!is.na(norm_evaluation_result$result),]

norm_wide = norm_evaluation_result %>% spread(norm_evaluation,result)

norm_sm = DE_methods(as.data.frame(norm_wide),res_column="silhouette_mean",met_column="norm_method")
impute_sm = DE_methods(as.data.frame(norm_wide),res_column="silhouette_mean",met_column="impute_method")

norm_var= aggregate(silhouette_mean ~ norm_method+cell_conditions, data=norm_wide, var)
norm_var= aggregate(silhouette_mean ~ norm_method, data=norm_var, mean)
norm_max= aggregate(silhouette_mean ~ norm_method+cell_conditions, data=norm_wide, max)
norm_max= aggregate(silhouette_mean ~ norm_method, data=norm_max, mean)


norm_sm_all = cbind(norm_max, norm_sm, norm_var,cluster_norm_sm, norm_trajectory_sm, norm_data_int_sm)
colnames(norm_sm_all) = c("norm_method", "sil_best","method","norm_coeff","p_value","norm_method1", "sil_var",
                          "method1","cluster_coeff","p_value",
                          "method2","trajectory_coeff","p_value",
                          "method3","data_int_coeff","p_value")
rownames(norm_sm_all) = norm_sm_all$method

norm_sm_all = norm_sm_all[,c("sil_best","norm_coeff","sil_var","cluster_coeff","trajectory_coeff","data_int_coeff")]
norm_sm_all = norm_sm_all[order(norm_sm_all$norm_coeff,decreasing = TRUE),]

##time
norm_timing_wide = norm_timing_result %>% spread(timing_method,result)
norm_timing_wide = norm_timing_wide[!is.na(norm_timing_wide$method_time),]
norm_timing_wide$method_time = norm_timing_wide$total_time-norm_timing_wide$method_time
norm_time_coeff = regress_running_time(as.data.frame(norm_timing_wide),"norm_method")

rownames(norm_time_coeff) = norm_time_coeff$method
norm_time_mean = aggregate(method_time ~ norm_method , data=norm_timing_wide, mean)
norm_time_mean$method_time = log2(norm_time_mean$method_time+1)
rownames(norm_time_mean) = norm_time_mean$norm_method 
time_df = cbind(norm_time_mean,norm_time_coeff)
time_df = time_df[,c("method_time","time_coeff")]
##


#### imputation

imp_var= aggregate(silhouette_mean ~ impute_method+cell_conditions, data=norm_wide, var)
imp_var= aggregate(silhouette_mean ~ impute_method, data=imp_var, mean)
imp_max= aggregate(silhouette_mean ~ impute_method+cell_conditions, data=norm_wide, max)
imp_max= aggregate(silhouette_mean ~ impute_method, data=imp_max, mean)

rownames(imp_max) = imp_max$impute_method
rownames(impute_sm) = impute_sm$method
rownames(imp_var) = imp_var$impute_method
rownames(cluster_impute_sm) = cluster_impute_sm$method
rownames(impute_trajectory_sm) = impute_trajectory_sm$method
rownames(impute_data_int_sm) = impute_data_int_sm$method

imp_sm_all = cbind(imp_max, impute_sm[rownames(imp_max),], 
                   imp_var[rownames(imp_max),],
                   cluster_impute_sm[rownames(imp_max),],
                   impute_trajectory_sm[rownames(imp_max),],
                   impute_data_int_sm[rownames(imp_max),])
colnames(imp_sm_all) = c("impute_method", "sil_best","method","norm_coeff","p_value","norm_method1", "sil_var",
                          "method1","cluster_coeff","p_value",
                          "method2","trajectory_coeff","p_value",
                          "method3","data_int_coeff","p_value")
rownames(imp_sm_all) = imp_sm_all$method
imp_sm_all = imp_sm_all[,c("sil_best","norm_coeff","sil_var","cluster_coeff","trajectory_coeff","data_int_coeff")]
imp_sm_all = imp_sm_all[2:4,]
imp_sm_all = imp_sm_all[order(imp_sm_all$sil_best,decreasing = TRUE),]

## time
norm_timing_wide = norm_timing_result %>% spread(timing_method,result)
norm_timing_wide = norm_timing_wide[!is.na(norm_timing_wide$method_time),]
#norm_timing_wide$method_time = norm_timing_wide$total_time-norm_timing_wide$method_time
norm_time_coeff = regress_running_time(as.data.frame(norm_timing_wide),"impute_method")

rownames(norm_time_coeff) = norm_time_coeff$method
norm_time_mean = aggregate(method_time ~ impute_method , data=norm_timing_wide, mean)
norm_time_mean$method_time = log2(norm_time_mean$method_time+1)
rownames(norm_time_mean) = norm_time_mean$impute_method 
imp_time_df = cbind(norm_time_mean,norm_time_coeff)
imp_time_df = imp_time_df[,c("method_time","time_coeff")]
##


norm_impute_df = rbind(scale(norm_sm_all),scale(imp_sm_all))
norm_impute_time = rbind(time_df, imp_time_df)
colnames(norm_impute_time) = c("Speed", "Scalability")
colnames(norm_impute_df) = c("Best performance", "Average performance", "Variability", 
                             "Performance in clustering",
                             "Performance in trajectory",
                             "Performance in data integration")

norm_impute_time$Speed = classify_speed(norm_impute_time$Speed)
norm_impute_time$Scalability = classify_scalability(norm_impute_time$Scalability)


pdf("summary_norm_imputation.pdf",width = 5,height = 6)
pheatmap(norm_impute_df,scale = "none",cluster_rows=FALSE,cluster_cols=FALSE,annotation_row=norm_impute_time,
         treeheight_row = 0, treeheight_col = 0,annotation_colors=anno_color,
         legend=FALSE,show_rownames=TRUE,show_colnames=TRUE,
         color=colorRampPalette(c(RColorBrewer::brewer.pal(3,"Spectral")[3], 
                                  RColorBrewer::brewer.pal(3,"Spectral")[2], 
                                  RColorBrewer::brewer.pal(3,"Spectral")[1]))(50),
         gaps_row=c(9),
         annotation_legend=FALSE,
         border_color="white",
         fontsize_row=10)
dev.off()




###### coefficient
DE_methods = function(res_spread,met_col = "traj_method"){
  res_method = c()
  res_m_type = c()
  res_coeff = c()
  res_pval = c()
  for(i in unique(as.character(res_spread$impute_method))){
    res_spread$the_method = "NO"
    res_spread$the_method[res_spread$impute_method==i] = "YES"
    fit = lm(result ~ the_method+cell_conditions,data=res_spread)
    sm = summary(fit)
    res_method = c(res_method,i)
    res_m_type = c(res_m_type,"impute_method")
    res_coeff = c(res_coeff,sm$coefficients[2,1])
    res_pval= c(res_pval,sm$coefficients[2,4])
  }
  for(i in unique(as.character(res_spread$norm_method))){
    res_spread$the_method = "NO"
    res_spread$the_method[res_spread$norm_method==i] = "YES"
    fit = lm(result ~ the_method+cell_conditions,data=res_spread)
    sm = summary(fit)
    res_method = c(res_method,i)
    res_m_type = c(res_m_type,"norm_method")
    res_coeff = c(res_coeff,sm$coefficients[2,1])
    res_pval= c(res_pval,sm$coefficients[2,4])
  }
  for(i in unique(as.character(res_spread[,met_col]))){
    res_spread$the_method = "NO"
    res_spread$the_method[res_spread[,met_col]==i] = "YES"
    fit = lm(result ~ the_method+cell_conditions,data=res_spread)
    sm = summary(fit)
    res_method = c(res_method,i)
    res_m_type = c(res_m_type,met_col)
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


DE_res = DE_methods(as.data.frame(trajectory_corr),"traj_method")
DE_res = DE_res[order(DE_res$coefficient),]
DE_res = rbind(DE_res[DE_res$method_type=="norm_method",],DE_res[DE_res$method_type=="impute_method",],DE_res[DE_res$method_type=="traj_method",])
DE_res$method_type[DE_res$method_type=="traj_method"] = "trajectory"
DE_res$method_type[DE_res$method_type=="norm_method"] = "normalization"
DE_res$method_type[DE_res$method_type=="impute_method"] = "imputation"
DE_res$method = factor(DE_res$method, levels = unique(DE_res$method))
pdf("coeff_corr_trajectory.pdf",width = 11,height = 5)
ggplot(data=DE_res,aes(x=method,y=coefficient,fill=method_type))+
  geom_bar(stat = "identity")+
  #scale_fill_brewer(palette="Set2")+
  scale_fill_manual(values = c("normalization"=brewer.pal(5, "Set2")[3],
                               "imputation"=brewer.pal(5, "Set2")[2],
                               "trajectory"=brewer.pal(5, "Set2")[5]))+
  ggtitle("Correlations")+
  labs(fill="")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

trajectory_over = trajectory_evaluation_result[trajectory_evaluation_result$traj_eval_method=="traj_overlap",]
trajectory_over = trajectory_over[!is.na(trajectory_over$result),]
trajectory_over$cell_conditions = trajectory_over$design

DE_res = DE_methods(as.data.frame(trajectory_over),"traj_method")
DE_res = DE_res[order(DE_res$coefficient),]
DE_res = rbind(DE_res[DE_res$method_type=="norm_method",],DE_res[DE_res$method_type=="impute_method",],DE_res[DE_res$method_type=="traj_method",])
DE_res$method_type[DE_res$method_type=="traj_method"] = "trajectory"
DE_res$method_type[DE_res$method_type=="norm_method"] = "normalization"
DE_res$method_type[DE_res$method_type=="impute_method"] = "imputation"
DE_res$method = factor(DE_res$method, levels = unique(DE_res$method))
pdf("coeff_over_trajectory.pdf",width = 11,height = 5)
ggplot(data=DE_res,aes(x=method,y=coefficient,fill=method_type))+
  geom_bar(stat = "identity")+
  #scale_fill_brewer(palette="Set2")+
  scale_fill_manual(values = c("normalization"=brewer.pal(5, "Set2")[3],
                               "imputation"=brewer.pal(5, "Set2")[2],
                               "trajectory"=brewer.pal(5, "Set2")[5]))+
  ggtitle("Overlaps")+
  labs(fill="")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

## data int

DE_res = DE_methods(as.data.frame(data_int_sil),"data_int_method")
DE_res = DE_res[order(DE_res$coefficient),]
DE_res = rbind(DE_res[DE_res$method_type=="norm_method",],DE_res[DE_res$method_type=="impute_method",],DE_res[DE_res$method_type=="data_int_method",])
DE_res$method_type[DE_res$method_type=="data_int_method"] = "data integration"
DE_res$method_type[DE_res$method_type=="norm_method"] = "normalization"
DE_res$method_type[DE_res$method_type=="impute_method"] = "imputation"
DE_res$method = factor(DE_res$method, levels = unique(DE_res$method))
pdf("coeff_sil_data_int.pdf",width = 11,height = 5)
ggplot(data=DE_res,aes(x=method,y=coefficient,fill=method_type))+
  geom_bar(stat = "identity")+
  #scale_fill_brewer(palette="Set2")+
  scale_fill_manual(values = c("normalization"=brewer.pal(5, "Set2")[3],
                               "imputation"=brewer.pal(5, "Set2")[2],
                               "data integration"=brewer.pal(5, "Set2")[4]))+
  ggtitle("Silhouette")+
  labs(fill="")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()


data_int_KBET = res_int_eval_all[res_int_eval_all$data_int_eval=="KBET",]
data_int_KBET = data_int_KBET[!(data_int_KBET$data_int_method=="zinbwave"),]
#data_int_sil = data_int_sil[!(data_int_sil$data_int_method=="no_int"),]
data_int_KBET = data_int_KBET[!is.na(data_int_KBET$result),]
data_int_KBET$cell_conditions = data_int_KBET$data


DE_res = DE_methods(as.data.frame(data_int_KBET),"data_int_method")
DE_res = DE_res[order(DE_res$coefficient),]
DE_res = rbind(DE_res[DE_res$method_type=="norm_method",],DE_res[DE_res$method_type=="impute_method",],DE_res[DE_res$method_type=="data_int_method",])
DE_res$method_type[DE_res$method_type=="data_int_method"] = "data integration"
DE_res$method_type[DE_res$method_type=="norm_method"] = "normalization"
DE_res$method_type[DE_res$method_type=="impute_method"] = "imputation"
DE_res$method = factor(DE_res$method, levels = unique(DE_res$method))
pdf("coeff_kBET_data_int.pdf",width = 11,height = 5)
ggplot(data=DE_res,aes(x=method,y=coefficient,fill=method_type))+
  geom_bar(stat = "identity")+
  #scale_fill_brewer(palette="Set2")+
  scale_fill_manual(values = c("normalization"=brewer.pal(5, "Set2")[3],
                                 "imputation"=brewer.pal(5, "Set2")[2],
                                 "data integration"=brewer.pal(5, "Set2")[4]))+
  ggtitle("kBET")+
  labs(fill="")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
