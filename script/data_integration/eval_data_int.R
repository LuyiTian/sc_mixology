#eval data integration
library(CellBench)
library(ggplot2)
library(scater)
library(scran)
library(kBET)
library(cluster)
library(plyr)
library(tidyr) 
library(ggrepel)
library(tidyverse)
library(bbplot)

silhouette_dr = function(sce){
  if("demuxlet_cls" %in% colnames(colData(sce))){
    sce = sce[,(sce$demuxlet_cls == "SNG")]
  }
  if("int_expr" %in% assayNames(sce)){
    try_res = try({
      sce = runPCA(sce,exprs_values="int_expr")
    })
    if (class(try_res) == "try-error") {
      return(NA)
    }
  }else{
    if(is.null(reducedDim(sce))){
      return(NA)
    }
  }
  ret_val = NA
  try_res = try({
    if("cell_line_demuxlet" %in% colnames(colData(sce))){
      sil = silhouette(as.numeric(factor(sce$cell_line_demuxlet)),  dist(reducedDim(sce)))
      ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))
    }else if("group" %in% colnames(colData(sce))){
      sil = silhouette(as.numeric(factor(sce$group)),  dist(reducedDim(sce)))
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


silhouette_dr_batch = function(sce){
  if("demuxlet_cls" %in% colnames(colData(sce))){
    sce = sce[,(sce$demuxlet_cls == "SNG")]
  }
  if("int_expr" %in% assayNames(sce)){
    try_res = try({
      sce = runPCA(sce,exprs_values="int_expr")
    })
    if (class(try_res) == "try-error") {
      return(NA)
    }
  }else{
    if(is.null(reducedDim(sce))){
      return(NA)
    }
  }
  ret_val = NA
  try_res = try({
    if("batch" %in% colnames(colData(sce))){
      sil = silhouette(as.numeric(factor(sce$batch)),  dist(reducedDim(sce)))
      ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))
    }
  })
  if (class(try_res) == "try-error") {
    return(NA)
  }
  return(ret_val)
}


KBET_dr = function(sce){
  if("cell_line_demuxlet" %in% colnames(colData(sce))){
    sce = sce[,(sce$cell_line_demuxlet %in% c("H1975","H2228","HCC827")) & (sce$demuxlet_cls == "SNG")]
  }
  if("int_expr" %in% assayNames(sce)){
    try_res = try({
      sce = runPCA(sce,exprs_values="int_expr",ncomponents = 2,feature_set=rownames(sce))
    })
    if (class(try_res) == "try-error") {
      return(NA)
    }
  }else{
    if(is.null(reducedDim(sce))){
      return(NA)
    }
  }
  ret_val = NA
  try_res = try({
    if("batch" %in% colnames(colData(sce))){
      batch.estimate <- kBET(t(reducedDim(sce)), batch=sce$batch,do.pca=FALSE, heuristic=TRUE, n_repeat=20,plot=FALSE)
      ret_val=1-batch.estimate$summary$kBET.observed[1]
    }
  })
  if (class(try_res) == "try-error") {
    return(NA)
  }
  return(ret_val)
}



data_int_eval = list(silhouette_dr=silhouette_dr,
                     silhouette_dr_batch=silhouette_dr_batch,
                     KBET=KBET_dr)


####
get_method_time = function(sce){
  typ="data_integration"
  if(!("int_expr" %in% assayNames(sce)) & length(reducedDimNames(sce))==0){
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
  typ="data_integration"
  if(!("int_expr" %in% assayNames(sce)) & length(reducedDimNames(sce))==0){
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
  for (i in unique(res$data_int_method)){
    tmp = res[res$data_int_method==i,]
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




RNAmix_after_data_int <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/RNAmix_after_data_int.Rds")

#res_int_eval = RNAmix_after_data_int %>% apply_methods(data_int_eval)
#saveRDS(res_int_eval,file="res_int_eval_RNAmix.Rds")

clu_time = RNAmix_after_data_int %>% apply_methods(timing_method)

for(i in 1:9){
  tmp = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/from_mac/RNAmix_after_scanorama_",i ,".Rds"))
  t_tmp = tmp %>% apply_methods(timing_method)
  clu_time = rbind(clu_time,t_tmp)
  #tmp = tmp %>% apply_methods(data_int_eval)
  #res_int_eval = rbind(res_int_eval,tmp)
}
for(i in 1:8){
  tmp = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/from_mac/sc_after_scanorama_",i ,".Rds"))
  t_tmp = tmp %>% apply_methods(timing_method)
  clu_time = rbind(clu_time,t_tmp)
  #tmp = tmp %>% apply_methods(data_int_eval)
  #res_int_eval = rbind(res_int_eval,tmp)
  tmp = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/data_int/sc_all_after_data_int_",i ,".Rds"))
  t_tmp = tmp %>% apply_methods(timing_method)
  clu_time = rbind(clu_time,t_tmp)
  #tmp = tmp %>% apply_methods(data_int_eval)
  #res_int_eval = rbind(res_int_eval,tmp)
}

saveRDS(clu_time,file="int_timing_result.Rds")

saveRDS(res_int_eval,file="res_int_eval_all.Rds")

# clu_time <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/int_timing_result.Rds")
# res_int_eval <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/res_int_eval_all.Rds")

ggplot(data=res_int_eval[res_int_eval$data_int_eval=="silhouette_dr",],aes(x=data_int_method,y=result,col=norm_method))+
  geom_boxplot()+
  bbc_style()


res_sel_top = function(res,topn=3){
  top_res = NULL
  for (i in unique(res$data_int_method)){
    tmp = res[res$data_int_method==i,]
    tmp = tmp[order(tmp$silhouette_dr,decreasing = TRUE)[1:min(topn,nrow(tmp))],]
    if(is.null(top_res)){
      top_res=tmp
    }else{
      top_res = rbind(tmp,top_res)
    }
  }
  return(top_res)
}


res_int_eval_RNAmix = res_int_eval[res_int_eval$data=="RNAmix",]
res_int_wide = res_int_eval_RNAmix %>% spread(data_int_eval,result)
res_int_wide_top = res_sel_top(res_int_wide,topn=2)
res_int_wide_top$norm_impute = paste(res_int_wide_top$norm_method,res_int_wide_top$impute_method,sep="_")
res_int_wide_top$data_int_method = revalue(res_int_wide_top$data_int_method, c("no_int"="none"))
pdf("Silhouette_kBET_RNAmix.pdf",width = 6,height = 7)
ggplot(data=res_int_wide_top,aes(x=silhouette_dr,y=KBET,col=data_int_method,shape=data_int_method,label=norm_impute))+
  geom_point(size=3)+
  geom_text_repel()+
  scale_colour_brewer(palette="Dark2")+
  scale_shape_manual(values=1:nlevels(res_int_wide_top$data_int_method))+
  labs(x="Silhouette of RNAmix",y="kBET (acceptance rate)",col="",shape="")+
  #bbc_style()
  theme_bw()+
  theme(text = element_text(size=15),legend.position="top")
dev.off()


res_int_eval_sc = res_int_eval[res_int_eval$data=="single_cell",]
res_int_wide = res_int_eval_sc %>% spread(data_int_eval,result)
res_int_wide_top = res_sel_top(res_int_wide,topn=2)
res_int_wide_top$norm_impute = paste(res_int_wide_top$norm_method,res_int_wide_top$impute_method,sep="_")
res_int_wide_top$data_int_method = revalue(res_int_wide_top$data_int_method, c("no_int"="none"))
pdf("Silhouette_kBET_single_cell.pdf",width = 6,height = 7)
ggplot(data=res_int_wide_top,aes(x=silhouette_dr,y=KBET,col=data_int_method,shape=data_int_method,label=norm_impute))+
  geom_point(size=3)+
  geom_text_repel()+
  scale_colour_brewer(palette="Dark2")+
  scale_shape_manual(values=1:nlevels(res_int_wide_top$data_int_method))+
  labs(x="Silhouette of single_cell",y="kBET (acceptance rate)",col="",shape="")+
  #bbc_style()
  theme_bw()+
  theme(text = element_text(size=15),legend.position="top")
dev.off()

#####
res_int_wide = res_int_eval_RNAmix %>% spread(data_int_eval,result)
res_int_wide$row_id = 1:nrow(res_int_wide)
res_int_wide$norm_impute = paste(res_int_wide$norm_method,res_int_wide$impute_method,sep="_")
res_int_wide_top = res_sel_top(res_int_wide,topn=2)
res_int_wide$norm_impute[!(res_int_wide$row_id %in% res_int_wide_top$row_id)] = NA
res_int_wide$whether_top = "NO"
res_int_wide$whether_top[(res_int_wide$row_id %in% res_int_wide_top$row_id)] = "YES"
res_int_wide = res_int_wide %>% gather(data_int_eval, result, silhouette_dr, silhouette_dr_batch)

res_int_wide$norm_impute[res_int_wide$data_int_eval=="silhouette_dr_batch"] = NA

res_int_wide$data_int_eval[res_int_wide$data_int_eval=="silhouette_dr"]="RNAmix"
res_int_wide$data_int_eval[res_int_wide$data_int_eval=="silhouette_dr_batch"]="batch"
res_int_wide$data_int_method = revalue(res_int_wide$data_int_method, c("no_int"="none"))

res_int_wide = res_int_wide[!(res_int_wide$data_int_method=="Seurat"), ]
res_int_wide = res_int_wide[!(res_int_wide$data_int_method=="zinbwave"), ]

pdf("RNAmix_data_int.pdf",width = 10,height = 5)
ggplot()+
  geom_boxplot(data=res_int_wide,aes(x=data_int_method,y=result,col=data_int_eval),outlier.colour = NA, position = position_dodge(width = 0.8))+
  geom_point(aes(x=res_int_wide$data_int_method[res_int_wide$whether_top=="YES"],y=res_int_wide$result[res_int_wide$whether_top=="YES"],col=res_int_wide$data_int_eval[res_int_wide$whether_top=="YES"]),position=position_jitterdodge(dodge.width=0.8,jitter.width=0.1),show.legend = FALSE)+
  #geom_text_repel(size=4,aes(x=res_int_wide$data_int_method[!is.na(res_int_wide$norm_impute)],y=res_int_wide$result[!is.na(res_int_wide$norm_impute)],col=res_int_wide$data_int_eval[!is.na(res_int_wide$norm_impute)],label=res_int_wide$norm_impute[!is.na(res_int_wide$norm_impute)]))+
  labs(y="Silhouette",col="type")+
  scale_colour_brewer(palette="Dark2")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

#####
res_int_wide = res_int_eval_sc %>% spread(data_int_eval,result)
res_int_wide$row_id = 1:nrow(res_int_wide)
res_int_wide$norm_impute = paste(res_int_wide$norm_method,res_int_wide$impute_method,sep="_")
res_int_wide_top = res_sel_top(res_int_wide,topn=2)
res_int_wide$norm_impute[!(res_int_wide$row_id %in% res_int_wide_top$row_id)] = NA
res_int_wide$whether_top = "NO"
res_int_wide$whether_top[(res_int_wide$row_id %in% res_int_wide_top$row_id)] = "YES"
res_int_wide = res_int_wide %>% gather(data_int_eval, result, silhouette_dr, silhouette_dr_batch)

res_int_wide$norm_impute[res_int_wide$data_int_eval=="silhouette_dr_batch"] = NA

res_int_wide$data_int_eval[res_int_wide$data_int_eval=="silhouette_dr"]="cell line"
res_int_wide$data_int_eval[res_int_wide$data_int_eval=="silhouette_dr_batch"]="batch"
res_int_wide$data_int_method = revalue(res_int_wide$data_int_method, c("no_int"="none"))
res_int_wide = res_int_wide[!(res_int_wide$data_int_method=="Seurat"), ]
res_int_wide = res_int_wide[!(res_int_wide$data_int_method=="zinbwave"), ]

pdf("SC_data_int.pdf",width = 10,height = 5)
ggplot()+
  geom_boxplot(data=res_int_wide,aes(x=data_int_method,y=result,col=data_int_eval),outlier.colour = NA, position = position_dodge(width = 0.8))+
  geom_point(aes(x=res_int_wide$data_int_method[res_int_wide$whether_top=="YES"],y=res_int_wide$result[res_int_wide$whether_top=="YES"],col=res_int_wide$data_int_eval[res_int_wide$whether_top=="YES"]),position=position_jitterdodge(dodge.width=0.8,jitter.width=0.1),show.legend = FALSE)+
  #geom_text_repel(size=4,aes(x=res_int_wide$data_int_method[!is.na(res_int_wide$norm_impute)],y=res_int_wide$result[!is.na(res_int_wide$norm_impute)],col=res_int_wide$data_int_eval[!is.na(res_int_wide$norm_impute)],label=res_int_wide$norm_impute[!is.na(res_int_wide$norm_impute)]))+
  labs(y="Silhouette",col="type")+
  scale_colour_brewer(palette="Dark2")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()


###### RNAmix all
all_data_int = as.character(unique(RNAmix_after_data_int$data_int_method))
all_norm = as.character(unique(RNAmix_after_data_int$norm_method))
all_impute = as.character(unique(RNAmix_after_data_int$impute_method))

for (it in all_data_int){
  for (nm in all_norm){
    for (imp in all_impute){
      tmp = RNAmix_after_data_int[RNAmix_after_data_int$data_int_method==it & RNAmix_after_data_int$norm_method==nm & RNAmix_after_data_int$impute_method==imp,]$result
      if(length(tmp)>0){
        sce = tmp[[1]]
        is_PCA = FALSE
        if(is.null(reducedDim(sce))){
          if (!("int_expr" %in% assayNames(sce))){
            next
          }
          sce = runPCA(sce,exprs_values="int_expr")
          is_PCA = TRUE
        }
        dr = reducedDim(sce)
        
        if(is_PCA){
          p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$group))+
            geom_point(size=4,alpha=0.7,show.legend = FALSE)+
            #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
            scale_shape_manual(values=c(0,1,2,3,4,5,8))+
            labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=2),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=2),"%)"),col="batch",shape="group")+
            theme_bw()+
            theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())
        }
        else{
          p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$group))+
            geom_point(size=4,alpha=0.7,show.legend = FALSE)+
            #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
            scale_shape_manual(values=c(0,1,2,3,4,5,8))+
            labs(x="Dim1",y="Dim2",col="batch",shape="group")+
            theme_bw()+
            theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())
        }
        pdf(file.path("data_int_plot",paste(nm,imp,it,"RNAmix.pdf",sep="_")))
        print(p)
        dev.off()
      }
    }
  }
}


####### RNAmix
sce = RNAmix_after_data_int[RNAmix_after_data_int$data_int_method=="no_int" & RNAmix_after_data_int$norm_method=="scran" & RNAmix_after_data_int$impute_method=="knn_smooth2",]$result[[1]]
col <- grDevices::rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
if(is.null(reducedDim(sce))){
  sce = runPCA(sce,exprs_values="int_expr")
}
dr = reducedDim(sce)
pdf("scran_knn_smooth2_noint_nolegend.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$group))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  scale_shape_manual(values=c(0,1,2,3,4,5,8))+
  labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=2),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=2),"%)"),col="batch",shape="group")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

#pdf("scran_knn_smooth2_noint.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$group))+
  geom_point(size=4,alpha=0.7,show.legend = T)+
  #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  scale_shape_manual(values=c(0,1,2,3,4,5,8))+
  labs(x="PC1",y="PC2",col="batch",shape="group")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
#dev.off()


sce = RNAmix_after_data_int[RNAmix_after_data_int$data_int_method=="Seurat" & RNAmix_after_data_int$norm_method=="Linnorm" & RNAmix_after_data_int$impute_method=="SAVER",]$result[[1]]
col <- grDevices::rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
if(is.null(reducedDim(sce))){
  sce = runPCA(sce,exprs_values="int_expr")
}
dr = reducedDim(sce)
pdf("Linnorm_SAVER_Seurat_nolegend.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$group))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  scale_shape_manual(values=c(0,1,2,3,4,5,8))+
  labs(x="PC1",y="PC2",col="batch",shape="group")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()


sce = RNAmix_after_data_int[RNAmix_after_data_int$data_int_method=="scMerge_us" & RNAmix_after_data_int$norm_method=="DESeq2" & RNAmix_after_data_int$impute_method=="SAVER",]$result[[1]]
col <- grDevices::rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
if(is.null(reducedDim(sce))){
  sce = runPCA(sce,exprs_values="int_expr")
}
dr = reducedDim(sce)
pdf("DESeq2_SAVER_scMerge_us_nolegend.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$group))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  scale_shape_manual(values=c(0,1,2,3,4,5,8))+
  labs(x="PC1",y="PC2",col="batch",shape="group")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

sce = RNAmix_after_data_int[RNAmix_after_data_int$data_int_method=="MNNs" & RNAmix_after_data_int$norm_method=="none" & RNAmix_after_data_int$impute_method=="knn_smooth2",]$result[[1]]
col <- grDevices::rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
if(is.null(reducedDim(sce))){
  sce = runPCA(sce,exprs_values="int_expr")
}
dr = reducedDim(sce)
pdf("none_knn_smooth2_MNNs_us_nolegend.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$group))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  scale_shape_manual(values=c(0,1,2,3,4,5,8))+
  labs(x="PC1",y="PC2",col="batch",shape="group")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

sce = RNAmix_after_data_int[RNAmix_after_data_int$data_int_method=="no_int" & RNAmix_after_data_int$norm_method=="Linnorm" & RNAmix_after_data_int$impute_method=="no_impute",]$result[[1]]
col <- grDevices::rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
if(is.null(reducedDim(sce))){
  sce = runPCA(sce,exprs_values="int_expr")
}
dr = reducedDim(sce)
pdf("Linnorm_no_impute_no_int_nolegend.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$group))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  scale_shape_manual(values=c(0,1,2,3,4,5,8))+
  labs(x="PC1",y="PC2",col="batch",shape="group")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

sce = RNAmix_after_data_int[RNAmix_after_data_int$data_int_method=="scMerge_s" & RNAmix_after_data_int$norm_method=="BASiCS" & RNAmix_after_data_int$impute_method=="SAVER",]$result[[1]]
col <- grDevices::rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
if(is.null(reducedDim(sce))){
  sce = runPCA(sce,exprs_values="int_expr")
}
dr = reducedDim(sce)
pdf("BASiCS_SAVER_scMerge_s_nolegend.pdf")
ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$group))+
  geom_point(size=4,alpha=0.7,show.legend = FALSE)+
  #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
  scale_shape_manual(values=c(0,1,2,3,4,5,8))+
  labs(x="PC1",y="PC2",col="batch",shape="group")+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()
################


###### single cell all

for(i in 1:8){
  tmp_after_int = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/data_int/sc_all_after_data_int_",i ,".Rds"))
  all_data_int = as.character(unique(tmp_after_int$data_int_method))
  all_norm = as.character(unique(tmp_after_int$norm_method))
  all_impute = as.character(unique(tmp_after_int$impute_method))
  for (it in all_data_int){
    for (nm in all_norm){
      for (imp in all_impute){
        tmp = tmp_after_int[tmp_after_int$data_int_method==it & tmp_after_int$norm_method==nm & tmp_after_int$impute_method==imp,]$result
        if(length(tmp)>0){
          sce = tmp[[1]]
          is_PCA = FALSE
          if(is.null(reducedDim(sce))){
            if (!("int_expr" %in% assayNames(sce))){
              next
            }
            sce = runPCA(sce,exprs_values="int_expr")
            is_PCA = TRUE
          }
          dr = reducedDim(sce)
          
          if(is_PCA){
            p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
              geom_point(size=3,alpha=0.5,show.legend = FALSE)+
              #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
              scale_shape_manual(values=c(16,0,1,18,2))+
              labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=2),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=2),"%)"),col="batch",shape="cell line")+
              theme_bw()+
              theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
          }
          else{
            p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
              geom_point(size=3,alpha=0.5,show.legend = FALSE)+
              #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
              scale_shape_manual(values=c(16,0,1,18,2))+
              labs(x="Dim1",y="Dim2",col="batch",shape="cell line")+
              theme_bw()+
              theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
          }
          pdf(file.path("data_int_plot","single_cell",paste(nm,imp,it,"sc_PCA.pdf",sep="_")))
          print(p)
          dev.off()
        }
      }
    }
  }
}
for(i in 1:8){
  tmp_after_int = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/from_mac/sc_after_scanorama_",i ,".Rds"))
  all_data_int = as.character(unique(tmp_after_int$data_int_method))
  all_norm = as.character(unique(tmp_after_int$norm_method))
  all_impute = as.character(unique(tmp_after_int$impute_method))
  for (it in all_data_int){
    for (nm in all_norm){
      for (imp in all_impute){
        tmp = tmp_after_int[tmp_after_int$data_int_method==it & tmp_after_int$norm_method==nm & tmp_after_int$impute_method==imp,]$result
        if(length(tmp)>0){
          sce = tmp[[1]]
          is_PCA = FALSE
          if(is.null(reducedDim(sce))){
            if (!("int_expr" %in% assayNames(sce))){
              next
            }
            sce = runPCA(sce,exprs_values="int_expr")
            is_PCA = TRUE
          }
          dr = reducedDim(sce)
          
          if(is_PCA){
            p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
              geom_point(size=3,alpha=0.5,show.legend = FALSE)+
              #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
              scale_shape_manual(values=c(16,0,1,18,2))+
              labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1]*100,digits=2),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2]*100,digits=2),"%)"),col="batch",shape="cell line")+
              theme_bw()+
              theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
          }
          else{
            p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
              geom_point(size=3,alpha=0.5,show.legend = FALSE)+
              #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
              scale_shape_manual(values=c(16,0,1,18,2))+
              labs(x="Dim1",y="Dim2",col="batch",shape="cell line")+
              theme_bw()+
              theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
          }
          pdf(file.path("data_int_plot","single_cell",paste(nm,imp,it,"sc_PCA.pdf",sep="_")))
          print(p)
          dev.off()
        }
      }
    }
  }
}

## single cell tSNE
for(i in 1:8){
  tmp_after_int = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/data_int/sc_all_after_data_int_",i ,".Rds"))
  all_data_int = as.character(unique(tmp_after_int$data_int_method))
  all_norm = as.character(unique(tmp_after_int$norm_method))
  all_impute = as.character(unique(tmp_after_int$impute_method))
  for (it in all_data_int){
    for (nm in all_norm){
      for (imp in all_impute){
        tmp = tmp_after_int[tmp_after_int$data_int_method==it & tmp_after_int$norm_method==nm & tmp_after_int$impute_method==imp,]$result
        if(length(tmp)>0){
          sce = tmp[[1]]
          #is_PCA = FALSE
          if(is.null(reducedDim(sce))){
            if (!("int_expr" %in% assayNames(sce))){
              next
            }
            sce = runTSNE(sce,exprs_values="int_expr",check_duplicates=FALSE)
            #is_PCA = TRUE
          }
          dr = reducedDim(sce)
          
            p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
              geom_point(size=3,alpha=0.5,show.legend = FALSE)+
              #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
              scale_shape_manual(values=c(16,0,1,18,2))+
              labs(x="Dim1",y="Dim2",col="batch",shape="cell line")+
              theme_bw()+
              theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
          pdf(file.path("data_int_plot","single_cell","TSNE",paste(nm,imp,it,"sc_TSNE.pdf",sep="_")))
          print(p)
          dev.off()
        }
      }
    }
  }
}
for(i in 1:8){
  tmp_after_int = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/from_mac/sc_after_scanorama_",i ,".Rds"))
  all_data_int = as.character(unique(tmp_after_int$data_int_method))
  all_norm = as.character(unique(tmp_after_int$norm_method))
  all_impute = as.character(unique(tmp_after_int$impute_method))
  for (it in all_data_int){
    for (nm in all_norm){
      for (imp in all_impute){
        tmp = tmp_after_int[tmp_after_int$data_int_method==it & tmp_after_int$norm_method==nm & tmp_after_int$impute_method==imp,]$result
        if(length(tmp)>0){
          sce = tmp[[1]]
          #is_PCA = FALSE
          if(is.null(reducedDim(sce))){
            if (!("int_expr" %in% assayNames(sce))){
              next
            }
            sce = runTSNE(sce,exprs_values="int_expr")
            #is_PCA = TRUE
          }
          dr = reducedDim(sce)
          
            p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
              geom_point(size=3,alpha=0.5,show.legend = FALSE)+
              #scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
              scale_shape_manual(values=c(16,0,1,18,2))+
              labs(x="Dim1",y="Dim2",col="batch",shape="cell line")+
              theme_bw()+
              theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
          pdf(file.path("data_int_plot","single_cell","TSNE",paste(nm,imp,it,"sc_TSNE.pdf",sep="_")))
          print(p)
          dev.off()
        }
      }
    }
  }
}

##### single cell

plot_df = res_int_wide_top[c(1,4,7,8,10,11,13,14),]
plot_df$norm_method = as.character(plot_df$norm_method)
plot_df$impute_method = as.character(plot_df$impute_method)
plot_df$data_int_method = as.character(plot_df$data_int_method)

set.seed(23333333)
for(i in 1:8){
  tmp = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/from_mac/sc_after_scanorama_",i ,".Rds"))
  tmp$norm_method = as.character(tmp$norm_method)
  tmp$impute_method = as.character(tmp$impute_method)
  tmp$data_int_method = as.character(tmp$data_int_method)
  for ( j in 1:nrow(plot_df)){
    tmp_sel = tmp[tmp$norm_method==plot_df$norm_method[j] & tmp$impute_method==plot_df$impute_method[j] & tmp$data_int_method==plot_df$data_int_method[j],]
    if (nrow(tmp_sel)>0){
      sce = tmp_sel$result[[1]]
      if(is.null(reducedDim(sce))){
        sce = runPCA(sce,exprs_values="int_expr")
      }
      dr = reducedDim(sce)
      if (plot_df$impute_method[j]=="Seurat"){
        lx = "Dim1"
        ly = "Dim2"
      }else{
        lx = "PC1"
        ly = "PC2"
      }
      pdf(paste(plot_df$data_int_method[j],plot_df$norm_method[j],plot_df$impute_method[j],"sc_int.pdf",sep="_"))
      p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
        geom_point(size=3,alpha=0.5,show.legend = FALSE)+
        labs(x=lx,y=ly,col="batch",shape="cell line")+
        scale_shape_manual(values=c(16,0,1,18,2))+
        theme_bw()+
        theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
      print(p)
      dev.off()

      sce = tmp_sel$result[[1]]
      if(is.null(reducedDim(sce))){
        sce = runTSNE(sce,exprs_values="int_expr")
        lx = "Dim1"
        ly = "Dim2"
      }else{
        if (plot_df$impute_method[j]=="Seurat"){
          lx = "Dim1"
          ly = "Dim2"
        }else{
          lx = "PC1"
          ly = "PC2"
        }
      }
      dr = reducedDim(sce)
      
      pdf(paste(plot_df$data_int_method[j],plot_df$norm_method[j],plot_df$impute_method[j],"sc_int_tsne.pdf",sep="_"))
      p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
        geom_point(size=3,alpha=0.5,show.legend = FALSE)+
        labs(x=lx,y=ly,col="batch",shape="cell line")+
        scale_shape_manual(values=c(16,0,1,18,2))+
        theme_bw()+
        theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
      print(p)
      dev.off()
    }
  }
}


for(i in 1:8){
  tmp = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/data_int/sc_all_after_data_int_",i ,".Rds"))
  tmp$norm_method = as.character(tmp$norm_method)
  tmp$impute_method = as.character(tmp$impute_method)
  tmp$data_int_method = as.character(tmp$data_int_method)
  for ( j in 1:nrow(plot_df)){
    tmp_sel = tmp[tmp$norm_method==plot_df$norm_method[j] & tmp$impute_method==plot_df$impute_method[j] & tmp$data_int_method==plot_df$data_int_method[j],]
    if (nrow(tmp_sel)>0){
      sce = tmp_sel$result[[1]]
      if(is.null(reducedDim(sce))){
        sce = runPCA(sce,exprs_values="int_expr")
      }
      dr = reducedDim(sce)
      if (plot_df$impute_method[j]=="Seurat"){
        lx = "Dim1"
        ly = "Dim2"
      }else{
        lx = "PC1"
        ly = "PC2"
      }
      pdf(paste(plot_df$data_int_method[j],plot_df$norm_method[j],plot_df$impute_method[j],"sc_int.pdf",sep="_"))
      p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
        geom_point(size=3,alpha=0.5,show.legend = FALSE)+
        labs(x=lx,y=ly,col="batch",shape="cell line")+
        scale_shape_manual(values=c(16,0,1,18,2))+
        theme_bw()+
        theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
      print(p)
      dev.off()
      
      sce = tmp_sel$result[[1]]
      if(is.null(reducedDim(sce))){
        sce = runTSNE(sce,exprs_values="int_expr")
        lx = "Dim1"
        ly = "Dim2"
      }else{
        if (plot_df$impute_method[j]=="Seurat"){
          lx = "Dim1"
          ly = "Dim2"
        }else{
          lx = "PC1"
          ly = "PC2"
        }
      }
      dr = reducedDim(sce)

      pdf(paste(plot_df$data_int_method[j],plot_df$norm_method[j],plot_df$impute_method[j],"sc_int_tsne.pdf",sep="_"))
      p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
        geom_point(size=3,alpha=0.5,show.legend = FALSE)+
        labs(x=lx,y=ly,col="batch",shape="cell line")+
        scale_shape_manual(values=c(16,0,1,18,2))+
        theme_bw()+
        theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
      print(p)
      dev.off()
    }
  }
}

pdf("sc_int_legend.pdf")
p = ggplot(data=NULL,aes(x=dr[,1],y=dr[,2],col=sce$batch,shape=sce$cell_line_demuxlet))+
  geom_point(size=3,alpha=0.5,show.legend = TRUE)+
  labs(x=lx,y=ly,col="batch",shape="cell line")+
  scale_shape_manual(values=c(16,0,1,18,2))+
  theme_bw()+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
print(p)
dev.off()



ggplot(data=as.data.frame(reducedDim(sce_sc_10x_5cl_qc,"PCA")),aes(x=PC1,y=PC3,shape=as.factor(sce_sc_10x_5cl_qc$demuxlet_cls),col=as.factor(colData(sce_sc_10x_5cl_qc)$cell_line_demuxlet)))+
  geom_point(size=3,show.legend = F,alpha=0.5)+
  scale_color_manual(values=c("H1975"="red","H2228" ="blue","HCC827"="green","H8383"="mediumorchid1", "A549"="yellow3"))+
  scale_shape_manual(values=c(0,16))+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())








sce = tmp$result[[3]]

sce = runPCA(sce,exprs_values="int_expr")
dr = reducedDim(sce)


sce = RNAmix_after_data_int[RNAmix_after_data_int$norm_method=="DESeq2" & RNAmix_after_data_int$impute_method=="SAVER" & RNAmix_after_data_int$data_int_method=="MNNs",]$result[[1]]
kbet_est <- KBET_dr(sce)
kbet_est

plotPCA(sce,run_args=list(exprs_values="int_expr"),colour_by="batch")
ggplot(data=NULL,aes(x=reducedDim(sce)[,1],y=reducedDim(sce)[,2],col=sce$group))+
  geom_point()+
  theme_bw()

batch.estimate <- kBET(t(assay(sce,"int_expr")[rowData(sce)$hi_var,]), k0=40, sce$batch,do.pca=TRUE,dim.pca=5, heuristic=FALSE, n_repeat=50,plot=TRUE)

tmp = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/from_mac/RNAmix_after_scanorama_",4 ,".Rds"))
sce = tmp$result[[3]]
plotPCA(sce,run_args=list(exprs_values="logcounts"),colour_by="batch")


