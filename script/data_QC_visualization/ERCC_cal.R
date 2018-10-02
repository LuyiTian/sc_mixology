library(scPipe)
library(scran)
library(scater)
library(ggplot2)
library(ggpmisc)
load("~/Dropbox/research/benchmark/rdata/singlecell_qc.RData")
load("~/Dropbox/research/benchmark/rdata/mRNAmix_qc.RData")

colData(sce2_qc)$ERCC_count = colSums(counts(sce2_qc[isSpike(sce2_qc),]))
colData(sce4_qc)$ERCC_count = colSums(counts(sce4_qc[isSpike(sce4_qc),]))
colData(sce8_qc)$ERCC_count = colSums(counts(sce8_qc[isSpike(sce8_qc),]))

dat = data.frame(experiment=c(rep("mRNAmix",ncol(sce2_qc)),
                              rep("single cell",ncol(sce4_qc)),
                              rep("mRNAmix_SORTseq",ncol(sce8_qc))),
                 ERCC_count=c(colData(sce2_qc)$ERCC_count,
                              colData(sce4_qc)$ERCC_count,
                              colData(sce8_qc)$ERCC_count),
                 count_per_cell=c((colData(sce2_qc)$total_count_per_cell),
                                  (colData(sce4_qc)$total_count_per_cell),
                                  (colData(sce8_qc)$total_count_per_cell)),
                 mRNA_amount=c(colData(sce2_qc)$mRNA_amount,
                               rep(NA,ncol(sce4_qc)),
                               colData(sce8_qc)$mRNA_amount))

dat = dat[dat$ERCC_count>2^(5.5),]


ggplot(data=dat[!(dat$experiment=="sincell_CELseq2"),],aes(x=experiment,y=ERCC_count,fill=factor(mRNA_amount)))+
  geom_boxplot()+
  theme_bw()


p1 = ggplot(data=dat[!(dat$experiment=="mRNAmix_SORTseq"),],aes(x=experiment,y=log2(ERCC_count),fill=factor(mRNA_amount)))+
  geom_boxplot()+
  scale_fill_brewer(palette="Spectral",direction = -1)+
  labs(fill="RNA amount (pg)",title="ERCC count distributions")+
  theme_bw()+
  theme(text = element_text(size=15,face="bold"),
        axis.text.x=element_text(angle = 30,hjust=1),
        axis.title.x=element_blank())

p12 =ggplot(data=dat[!(dat$experiment=="mRNAmix_SORTseq"),],aes(x=experiment,y=log2(count_per_cell),fill=factor(mRNA_amount)))+
  geom_boxplot()+
  scale_fill_brewer(palette="Spectral",direction = -1)+
  labs(fill="RNA amount (pg)",title="total count distributions")+
  theme_bw()+
  theme(text = element_text(size=15,face="bold"),
        axis.text.x=element_text(angle = 30,hjust=1),
        axis.title.x=element_blank())


p2 = ggplot(data=dat[!(dat$experiment=="mRNAmix_SORTseq"),],aes(x=log2(count_per_cell),y=log2(ERCC_count),col=experiment))+
  geom_point(alpha=0.9)+geom_smooth(method='lm')+
  labs(title="count relationships in both")+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,sep = "~~~~")), 
               label.x.npc = "left", label.y.npc = 0.95, formula = y ~ x  ,parse = TRUE, size = 3)+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  theme(text = element_text(size=15,face="bold"))





ggplot(data=dat[,],aes(x=count_per_cell,y=log2(ERCC_count),col=experiment))+
  geom_point(alpha=0.7)+geom_smooth(method='lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,sep = "~~~~")), 
               label.x.npc = "right", label.y.npc = 0.05, formula = y ~ x  ,parse = TRUE, size = 3)+
  theme_bw()



p3=ggplot(data=dat[(dat$experiment=="mRNAmix"),],aes(x=log2(count_per_cell),y=log2(ERCC_count),col=factor(mRNA_amount)))+
  geom_point(alpha=0.9)+geom_smooth(method='lm')+
  scale_color_brewer(palette="Spectral",direction = -1)+
  labs(title="count relationships in mRNAmix",col="RNA amount (pg)")+
  stat_poly_eq(aes(label = paste( ..adj.rr.label..)), 
               label.x.npc = "left", show.legend=T,label.y.npc = 0.95, formula = y ~ x  ,parse = TRUE, size = 3)+
  theme_bw()+
  theme(text = element_text(size=15,face="bold"))


ggplot(data=dat[(dat$experiment=="mRNAmix_SORTseq"),],aes(x=count_per_cell,y=log2(ERCC_count),col=factor(mRNA_amount)))+
  geom_point(alpha=0.7)+geom_smooth(method='lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,sep = "~~~~")), 
               label.x.npc = "right", label.y.npc = 0.15, formula = y ~ x  ,parse = TRUE, size = 3)+
  theme_bw()





#####
tmpp = dat[(dat$experiment=="mRNAmix"),]
tmpp$ERCC_count = log2(tmpp$ERCC_count)
tmpp$count_per_cell = log2(tmpp$count_per_cell)

pdf("ERCC_count_mRNAmix.pdf")
marginal_plot(x = count_per_cell, y = ERCC_count, group = mRNA_amount, 
              data = tmpp, bw = "nrd", lm_formula =  y ~ x , pch = 16, cex = 1,lm_show = TRUE)
dev.off()

tmpp = dat[(dat$experiment=="mRNAmix"),]
tmpp$ERCC_count = log2(tmpp$ERCC_count)
tmpp$count_per_cell = log2(tmpp$count_per_cell)

pdf("ERCC_count_mRNAmix_nogroup.pdf")
marginal_plot(x = count_per_cell, y = ERCC_count, 
              data = tmpp, bw = "nrd", lm_formula =  y ~ x , pch = 16, cex = 1,lm_show = TRUE)
dev.off()

tmpp = dat[!(dat$experiment=="mRNAmix_SORTseq"),]
tmpp$ERCC_count = log2(tmpp$ERCC_count)
tmpp$count_per_cell = log2(tmpp$count_per_cell)

pdf("ERCC_count_mRNAmix_sincell.pdf")
marginal_plot(x = count_per_cell, y = ERCC_count, group = experiment, 
              data = tmpp, bw = "nrd", lm_formula =  y ~ x , pch = 16, cex = 1,lm_show = TRUE)
dev.off()

#####


####

library(edgeR)
bulk.count.matrix <- read.delim("~/Dropbox/research/benchmark/reverse.stranded.unfiltered.count.matrix.txt", row.names=1, stringsAsFactors=FALSE)
bulk.count.matrix = bulk.count.matrix[,c("RSCE_6_BC2CTUACXX_CTTGTA_L005_R1.bam", "RSCE_8_BC2CTUACXX_CCGTCC_L007_R1.bam",
                                         "RSCE_10_BC2CTUACXX_AGTTCC_L006_R1.bam", "RSCE_12_BC2CTUACXX_TGACCA_L005_R1.bam",
                                         "RSCE_14_BC2CTUACXX_GTCCGC_L003_R1.bam", "RSCE_16_BC2CTUACXX_AGTTCC_L004_R1.bam")]
colnames(bulk.count.matrix) = c("H1975_rep1", "H1975_rep2",
                                "HCC827_rep1", "HCC827_rep2",
                                "H2228_rep1", "H2228_rep2")
grp = c("H1975","H1975","HCC827","HCC827","H2228","H2228")
y = DGEList(counts=bulk.count.matrix, group = grp)
y = estimateCommonDisp(y)

design <- model.matrix(~0+grp)
fit = glmFit(y,design=design)
contrasts = makeContrasts(
  H2228=grpH2228-0.5*grpH1975-0.5*grpHCC827,
  H1975=grpH1975-0.5*grpHCC827-0.5*grpH2228,
  HCC827=grpHCC827-0.5*grpH2228-0.5*grpH1975,
  levels=colnames(design))

lrt = glmLRT(fit, contrast=contrasts[,1])
H2228 = topTags(lrt,n=500,p.value=0.001)

lrt = glmLRT(fit, contrast=contrasts[,2])
H1975 = topTags(lrt,n=500,p.value=0.001)

lrt = glmLRT(fit, contrast=contrasts[,3])
HCC827 = topTags(lrt,n=500,p.value=0.001)




sce10x_qc = convert_geneid(sce10x_qc, returns = "entrezgene")



bulk.count.matrix = bulk.count.matrix[rownames(bulk.count.matrix) %in% c(rownames(H2228), rownames(H1975), rownames(HCC827)),]
uni_genes = intersect(rownames(bulk.count.matrix), rownames(sce10x_qc))


bulk_mat_sel = bulk.count.matrix[uni_genes,]


sce10x_qc = computeSumFactors(sce10x_qc)
sce10x_qc <- normalize(sce10x_qc)
sce_qc_sel = sce10x_qc[uni_genes,]
# H1975
H1975_corr = apply(logcounts(sce_qc_sel), 2, function(x){cor(bulk_mat_sel[,"H1975_rep1"], x, method="spearman")})
H2228_corr = apply(logcounts(sce_qc_sel), 2, function(x){cor(bulk_mat_sel[,"HCC827_rep1"], x, method="spearman")})
HCC827_corr = apply(logcounts(sce_qc_sel), 2, function(x){cor(bulk_mat_sel[,"H2228_rep1"], x, method="spearman")})
cor_mat = cbind(H1975_corr, H2228_corr, HCC827_corr)

library(mclust)
gmm_out = Mclust(cor_mat,3,modelNames="VVV",verbose = FALSE)
plot(gmm_out, what = "classification")


colData(sce10x_qc)$cla =gmm_out$classification

ggplot(data=as.data.frame(colData(sce4_qc)),aes(x=log2(total_count_per_cell),y=log2(ERCC_count),col=factor(cla)))+
  geom_point(alpha=0.7)+geom_smooth(method='lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,sep = "~~~~")), 
               label.x.npc = "right", label.y.npc = 0.15, formula = y ~ x  ,parse = TRUE, size = 3)+
  theme_bw()

ggplot(data=as.data.frame(colData(sce4_qc)),aes(x=log2(total_count_per_cell),y=log2(colData(sce4_qc)$ERCC_count),col=factor(cla)))+
  geom_point()+
  theme_bw()

ggplot(data=dat[(dat$experiment=="sincell_CELseq2"),],aes(x=count_per_cell,y=ERCC_count,col=factor(gmm_out$classification)))+
  geom_point(alpha=0.7)+geom_smooth(method='lm')+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,sep = "~~~~")), 
               label.x.npc = "right", label.y.npc = 0.15, formula = y ~ x  ,parse = TRUE, size = 3)+
  theme_bw()



dat$ERCC_count = log2(dat$ERCC_count)
dat$mRNA_amount = dat$mRNA_amount/3.75

aaa = as.data.frame(colData(sce2_qc)[,c("H2228_prop", "H1975_prop", "HCC827_prop", "ERCC_count", "mRNA_amount", "total_count_per_cell")])
aaa$total_count_per_cell = log2(aaa$total_count_per_cell)
aaa$ERCC_count = log2(aaa$ERCC_count)

fit = aov(formula = total_count_per_cell ~ H2228_prop +HCC827_prop+mRNA_amount+ERCC_count, 
         data = aaa)

summary(fit)
su = summary(fit)
sum_sq = su[[1]]$`Sum Sq`

x = c("mixture design", "mRNA amount", "ERCC counts", "Residuals")

p4 = ggplot(data=NULL,aes(x=x,y=c(sum(sum_sq[1:3]),sum_sq[4],sum_sq[5],sum_sq[6]),fill=x))+
  geom_bar(stat="identity",show.legend = F)+
  labs(fill="method",x="top", y="Sums of Squares")+
  theme_bw()+
  labs(title="Decomposing count variance")+
  theme(text = element_text(size=15,face="bold"),
        axis.text.x=element_text(angle = 30,hjust=1),
        axis.title.x=element_blank())



library(ggpubr)
pdf("ERCC_mRNA_sincell.pdf",width = 12,height = 8)
#ggarrange(p1,p2,p3,ncol = 1, nrow = 3,legend="right",labels="AUTO",font.label=list(size = 20, face = "bold"),common.legend=T)
ggarrange(p1,p4,p3,p2,ncol = 2, nrow = 2,labels="AUTO")
dev.off()


pdf("ERCC_mRNA_boxplot.pdf",width = 11,height = 4)
#ggarrange(p1,p2,p3,ncol = 1, nrow = 3,legend="right",labels="AUTO",font.label=list(size = 20, face = "bold"),common.legend=T)
ggarrange(p1,p12,ncol = 2, nrow = 1)
dev.off()

pdf("ERCC_mRNA_lm.pdf",width = 11,height = 4)
#ggarrange(p1,p2,p3,ncol = 1, nrow = 3,legend="right",labels="AUTO",font.label=list(size = 20, face = "bold"),common.legend=T)
ggarrange(p3,p2,ncol = 2, nrow = 1)
dev.off()


#############


colData(sce_9cells_qc)$ERCC_count = colSums(counts(sce_9cells_qc[isSpike(sce_9cells_qc),]))


sce_9cells_qc = detect_outlier(sce_9cells_qc,comp=2)

sce_9cells_rm = remove_outliers(sce_9cells_qc)

ggplot(data=as.data.frame(colData(sce_9cells_rm)),aes(y=log2(ERCC_count),
                                                      x=log2(total_count_per_cell),
                                                      col=factor(rowSums(as.data.frame(colData(sce_9cells_rm))[,c("H1975","H2228","HCC827")]))))+geom_point()

ggplot(data=as.data.frame(colData(sce_9cells_rm)),aes(y=log2(ERCC_count),
                                                      x=log2(total_count_per_cell),
                                                      col=poor_quality))+geom_point()

