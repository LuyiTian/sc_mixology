# plot QC stat

library(ggplot2)
library(scPipe)
library(ggpubr)
setwd("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit")

get_mapping_rate = function(sce){
  if("mapped_to_ERCC" %in% colnames(colData(sce))){
    mapping_rate = (sce$mapped_to_exon+sce$mapped_to_ERCC)/(sce$mapped_to_exon+sce$mapped_to_ERCC+sce$unaligned+sce$aligned_unmapped+sce$mapped_to_intron+sce$ambiguous_mapping)
  }else{
    mapping_rate = (sce$mapped_to_exon)/(sce$mapped_to_exon+sce$unaligned+sce$aligned_unmapped+sce$mapped_to_intron+sce$ambiguous_mapping)
  }
  return(mapping_rate)
}

get_mapping_reads = function(sce){
  if("mapped_to_ERCC" %in% colnames(colData(sce))){
    mapping_rate = (sce$mapped_to_exon+sce$mapped_to_ERCC)
  }else{
    mapping_rate = (sce$mapped_to_exon)
  }
  return(mapping_rate)
}


get_intron_mapping_rate = function(sce){
  if("mapped_to_ERCC" %in% colnames(colData(sce))){
    mapping_rate = (sce$mapped_to_intron)/(sce$mapped_to_exon+sce$mapped_to_ERCC+sce$unaligned+sce$aligned_unmapped+sce$mapped_to_intron+sce$ambiguous_mapping)
  }else{
    mapping_rate = (sce$mapped_to_intron)/(sce$mapped_to_exon+sce$unaligned+sce$aligned_unmapped+sce$mapped_to_intron+sce$ambiguous_mapping)
  }
  return(mapping_rate)
}

get_amp_rate = function(sce){
  if("mapped_to_ERCC" %in% colnames(colData(sce))){
    amp_rate = (sce$mapped_to_exon+sce$mapped_to_ERCC)/(sce$total_count_per_cell)
  }else{
    amp_rate = (sce$mapped_to_exon)/(sce$total_count_per_cell)
  }
  return(amp_rate)
}


load("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/sincell_with_class_5cl.RData")
load("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/sincell_with_class_new.RData")
load("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/CellBench_data/data/mRNAmix_qc.RData")
load("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/CellBench_data/data/9cellmix_qc.RData")

sce_l = list(sce_SC1_qc,
             sce_SC2_qc,
             sce_SC3_qc,
             sce_SC4_qc,
             sce_POP_sel_qc,
             sce2_qc,
             sce8_qc,
             sce_sc_10x_qc,
             sce_sc_CELseq2_qc,
             sce_sc_Dropseq_qc,
             sce_sc_10x_5cl_qc,
             sc_Celseq2_5cl_p1,
             sc_Celseq2_5cl_p2,
             sc_Celseq2_5cl_p3)

data_info = c(rep("cellmix1",ncol(sce_SC1_qc)),
              rep("cellmix2",ncol(sce_SC2_qc)),
              rep("cellmix3",ncol(sce_SC3_qc)),
              rep("cellmix4",ncol(sce_SC4_qc)),
              rep("cellmix5",ncol(sce_POP_sel_qc)),
              rep("RNAmix_CELseq2",ncol(sce2_qc)),
              rep("RNAmix_Sortseq",ncol(sce8_qc)),
              rep("sc_10X",ncol(sce_sc_10x_qc)),
              rep("sc_CELseq2",ncol(sce_sc_CELseq2_qc)),
              rep("sc_Dropseq",ncol(sce_sc_Dropseq_qc)),
              rep("sc_10x_5cl",ncol(sce_sc_10x_5cl_qc)),
              rep("sc_CELseq2_5cl_p1",ncol(sc_Celseq2_5cl_p1)),
              rep("sc_CELseq2_5cl_p2",ncol(sc_Celseq2_5cl_p2)),
              rep("sc_CELseq2_5cl_p3",ncol(sc_Celseq2_5cl_p3))
              )

mapping_rate = unlist(lapply(sce_l,get_mapping_rate))
df = data.frame(batch=data_info,mapping_rate=mapping_rate)
df$batch = factor(df$batch,levels = unique(data_info))

p1 = ggplot(data=df,aes(x=batch,y=mapping_rate,col=batch))+
  geom_boxplot(show.legend = FALSE)+
  labs(x="",y="exon mapping rate")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
p1

mapping_rate = unlist(lapply(sce_l,get_intron_mapping_rate))
df = data.frame(batch=data_info,mapping_rate=mapping_rate)
df$batch = factor(df$batch,levels = unique(data_info))

p2 = ggplot(data=df,aes(x=batch,y=mapping_rate,col=batch))+
  geom_boxplot(show.legend = FALSE)+
  labs(x="",y="intron mapping rate")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
p2


mapping_rate = unlist(lapply(sce_l,get_amp_rate))
df = data.frame(batch=data_info,mapping_rate=mapping_rate)
df$batch = factor(df$batch,levels = unique(data_info))

p3 = ggplot(data=df,aes(x=batch,y=log2(mapping_rate+1),col=batch))+
  geom_boxplot(show.legend = FALSE)+
  labs(x="",y="log2(amp_rate)")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
p3

mapping_rate = unlist(lapply(sce_l,get_mapping_reads))
df = data.frame(batch=data_info,mapping_rate=mapping_rate)
df$batch = factor(df$batch,levels = unique(data_info))

p4 = ggplot(data=df,aes(x=batch,y=log2(mapping_rate+1),col=batch))+
  geom_boxplot(show.legend = FALSE)+
  labs(x="",y="log2(exon_mapped_reads)")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
p4


mapping_rate = unlist(lapply(sce_l,function(x){x$total_count_per_cell}))
df = data.frame(batch=data_info,mapping_rate=mapping_rate)
df$batch = factor(df$batch,levels = unique(data_info))

p5 = ggplot(data=df,aes(x=batch,y=log2(mapping_rate+1),col=batch))+
  geom_boxplot(show.legend = FALSE)+
  labs(x="",y="log2(UMI_count_per_cell)")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
p5


mapping_rate = unlist(lapply(sce_l,function(x){x$number_of_genes}))
df = data.frame(batch=data_info,mapping_rate=mapping_rate)
df$batch = factor(df$batch,levels = unique(data_info))

p6 = ggplot(data=df,aes(x=batch,y=mapping_rate,col=batch))+
  geom_boxplot(show.legend = FALSE)+
  labs(x="",y="number_of_genes")+
  theme_bw()+
  theme(text = element_text(size=20),axis.title.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
p6



pdf("QC_summary_all.pdf",width = 14,height = 20)
ggarrange(p2,p1,p4,p5,p3,ncol = 1,nrow = 5,align="hv")
dev.off()


df_comp = data.frame(gene_number = unlist(lapply(sce_l,function(x){x$number_of_genes})), 
                     UMI_count = unlist(lapply(sce_l,function(x){x$total_count_per_cell})),
                     exon_count = unlist(lapply(sce_l,get_mapping_reads)),
                     batch = data_info,stringsAsFactors = FALSE)

df_comp = df_comp[grepl("sc_",df_comp$batch),]
df_comp$protocol = "10x"
df_comp$protocol[grepl("CEL",df_comp$batch)] = "CEL-seq2"
df_comp$protocol[grepl("Drop",df_comp$batch)] = "Drop-seq"

pdf("QC_metrics_protocols_benchmark_data.pdf")
ggplot(data=df_comp,aes(x=exon_count,y=gene_number,col=protocol))+geom_point(alpha=0.5)+theme_bw()+theme(text = element_text(size=15))

ggplot(data=df_comp,aes(x=UMI_count,y=gene_number,col=protocol))+geom_point(alpha=0.5)+theme_bw()+theme(text = element_text(size=15))

ggplot(data=df_comp,aes(x=exon_count,y=UMI_count,col=protocol))+geom_point(alpha=0.5)+theme_bw()+theme(text = element_text(size=15))

dev.off()






