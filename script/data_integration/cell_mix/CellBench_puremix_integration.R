## ----global_options, include=FALSE---------------------------------------
# key libraries
library(knitr)

# set options to knitting
eval_code = FALSE
to_run = FALSE

# global options

knitr::opts_chunk$set(dpi = 100, print = TRUE, echo=FALSE, warning=FALSE, message=FALSE, eval = TRUE, fig.show=TRUE, fig.width= 7,fig.height= 6,fig.align='center', out.width = '60%', fig.path= 'Figures-puremix/')

library(tictoc) # to record timing
library(kableExtra) # for a nice table

## ---- eval  = FALSE------------------------------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite('SingleCellExperiment')
## biocLite('scran')
## biocLite('scater')
## biocLite('Seurat')
## biocLite('zinbwave')
## 
## install.packages("BiocManager")
## BiocManager::install("zinbwave")
## 
## # for scMerge:
## # Some CRAN packages required by scMerge
## install.packages(c("ruv", "rsvd", "igraph", "pdist", "proxy", "foreach", "doSNOW", "distr"))
## devtools::install_github("theislab/kBET")
## 
## # Some BioConductor packages required by scMerge
## # try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R")
## biocLite(c("SingleCellExperiment", "M3Drop"))
## # Installing scMerge using
## devtools::install_github("SydneyBioX/scMerge")

## ------------------------------------------------------------------------
library(tictoc) # to record timing

library(SingleCellExperiment) 
# for normalisation:
library(scran)  # also for MNNcorrect
library(scater)
# for PCA 
library(mixOmics)  # by default centers but does not scale the data
# for seura CCA and seurat normalisation
library(Seurat)
library(zinbwave)
library(scRNAseq)  # to extract weights in zinbwave

library(scMerge)

library(reticulate)
scanorama <- import('scanorama')

library(kBET)

library(cluster) # to calculate silhouette and ARI
library(clues)   # to calculate ARI

## ---- include = FALSE----------------------------------------------------
# to run code already saved
if(!to_run & file.exists('CellBench_puremix.results.RData')) load('CellBench_puremix.results.RData')

# if to run then run the whole code
if(to_run) eval_code = TRUE

## ---- eval = eval_code---------------------------------------------------
## load('Data/sincell_with_class.RData')

## ------------------------------------------------------------------------
dim(counts(sce10x_qc))  # isolation type, QC'ed with scPipe, cell type class indicated
dim(counts(sce4_qc))
dim(counts(scedrop_qc_qc))

## ------------------------------------------------------------------------
dt = data.frame(
  Chromium10X = summary(as.factor(sce10x_qc$cell_line)),
  CELseq2 = summary(as.factor(sce4_qc$cell_line)),
  DROPseq = summary(as.factor(scedrop_qc_qc$cell_line))
)

kable(t(dt), caption = 'Number of cells per cell line type and per protocol')

## ---- eval = eval_code---------------------------------------------------
## sc10x.norm = computeSumFactors(sce10x_qc)
## sc10x.norm = normalize(sc10x.norm)
## 
## scdrop.norm = computeSumFactors(scedrop_qc_qc)
## scdrop.norm = normalize(scdrop.norm)
## # to deal with identical names during the sequencing, we rename the DROP-seq samples
## colnames(scdrop.norm) = paste0("dropseq_",colnames(scdrop.norm))
## 
## sccel.norm = computeSumFactors(sce4_qc)
## sccel.norm = normalize(sccel.norm)

## ---- eval = eval_code---------------------------------------------------
## # using mixOmics for PCA, by default centers but does not scale
## pca.res.10x = mixOmics::pca(t(logcounts(sc10x.norm)), ncomp = 3)
## pca.res.celseq = mixOmics::pca(t(logcounts(sccel.norm)), ncomp = 3)
## pca.res.dropseq = mixOmics::pca(t(logcounts(scdrop.norm)), ncomp = 3)
## 
## # the plot option to look at the explained variance per component:
## # plot(pca.res.10x)

## ----PCA-10X-------------------------------------------------------------
plotIndiv(pca.res.10x, pch = 1, 
          group = sce10x_qc$cell_line, col.per.group = color.mixo(4:6),
          legend = TRUE, title = 'PCA 10X' 
          )

## ----PCA-CELseq2---------------------------------------------------------
plotIndiv(pca.res.celseq, pch = 2, 
          group = sce4_qc$cell_line, col.per.group = color.mixo(4:6),
          legend = TRUE, title = 'PCA CEL-seq2')

## ----PCA-DROPseq---------------------------------------------------------
plotIndiv(pca.res.dropseq, pch = 3, 
          group = scedrop_qc_qc$cell_line, col.per.group = color.mixo(4:6),
          legend = TRUE, title = 'PCA DROP-seq')

## ------------------------------------------------------------------------
# intersection of the UMI
list.intersect = Reduce(intersect, list(rownames(logcounts(sc10x.norm)), rownames(logcounts(sccel.norm)), rownames(logcounts(scdrop.norm))))
#length(list.intersect)

## ---- eval = eval_code---------------------------------------------------
## data.combined = t(data.frame(logcounts(sc10x.norm)[list.intersect,], logcounts(sccel.norm)[list.intersect,], logcounts(scdrop.norm)[list.intersect,]))
## dim(data.combined)
## 
## # cell line information to assign to the combined data for plotting
## cell.line = as.factor(c(sce10x_qc$cell_line, sce4_qc$cell_line, scedrop_qc_qc$cell_line))
## names(cell.line) = rownames(data.combined)
## 
## # batch information (protocol)
## batch = as.factor(c(rep('10X', ncol(logcounts(sc10x.norm))), rep('CEL-seq2', ncol(logcounts(sccel.norm))),
##                   rep('DROP-seq', ncol(logcounts(scdrop.norm)))
##                     ))
## names(batch) = rownames(data.combined)

## ------------------------------------------------------------------------
kable(summary(batch), caption = 'Number of cells per protocol')

## ---- eval = eval_code---------------------------------------------------
## pca.combined = mixOmics::pca(data.combined, ncomp = 2)

## ----PCA-combined--------------------------------------------------------
# color indicates protocol
plotIndiv(pca.combined, pch = cell.line, group = batch, legend = TRUE, legend.title = 'Protocol', legend.title.pch = 'Cell line', title = 'PCA')

## ----Seurat-normalisation, eval = eval_code------------------------------
## srt10x <- CreateSeuratObject(raw.data = counts(sce10x_qc))
## srt10x <- NormalizeData(object = srt10x)
## srt10x <- ScaleData(object = srt10x)
## srt10x <- FindVariableGenes(object = srt10x, do.plot = FALSE, topn=10000, display.progress = FALSE)
## 
## srtCEL <- CreateSeuratObject(raw.data = counts(sce4_qc))
## srtCEL <- NormalizeData(object = srtCEL)
## srtCEL <- ScaleData(object = srtCEL)
## srtCEL <- FindVariableGenes(object = srtCEL, do.plot = FALSE, topn=10000, display.progress = FALSE)
## 
## colnames(scedrop_qc_qc) = paste0("dropseq_",colnames(scedrop_qc_qc))
## srtDROP <- CreateSeuratObject(raw.data = counts(scedrop_qc_qc))
## srtDROP <- NormalizeData(object = srtDROP)
## srtDROP <- ScaleData(object = srtDROP)
## srtDROP <- FindVariableGenes(object = srtDROP, do.plot = FALSE, topn=10000, display.progress = FALSE)
## 
## # extract meta data indicating the protocol for diagonal CCA
## srt10x@meta.data[, "protocol"] <- "10X"
## srt10x@meta.data[, "cell line"] <- sce10x_qc$cell_line
## srtCEL@meta.data[, "protocol"] <- "CELSeq"
## srtCEL@meta.data[, "cell line"] <- sce4_qc$cell_line
## srtDROP@meta.data[, "protocol"] <- "Drop-seq"
## srtDROP@meta.data[, "cell line"] <- scedrop_qc_qc$cell_line
## 
## # then calculate intersection between those lists
## high.var.genes = Reduce(intersect, list(srt10x@var.genes, srtCEL@var.genes, srtDROP@var.genes))
## length(high.var.genes)

## ---- eval = eval_code---------------------------------------------------
## scran_high_var = function(
##     sce,    # sce object
##     topn=2000  # by default we look at the top 3,000 genes
##     ){
##   var.fit <- trendVar(sce, method="loess", use.spikes=FALSE)
##   # calclates biological and technical variance
##   var.out <- decomposeVar(sce, var.fit)
##   # order genes with high biological variance.
##   hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
##   return(invisible(rownames(hvg.out)))
## }
## 
## sce10x.high.var = scran_high_var(sc10x.norm)
## sccel.high.var = scran_high_var(sccel.norm)
## scedrop.high.var = scran_high_var(scdrop.norm)
## 
## # then calculate intersection between those lists
## high.var.genes.scran = Reduce(intersect, list(sce10x.high.var, sccel.high.var, scedrop.high.var))
## length(high.var.genes.scran)

## ---- eval = eval_code---------------------------------------------------
## tic('MNN high variable genes')
## MNNcorrect.res.high.var = mnnCorrect(logcounts(sc10x.norm)[high.var.genes.scran,], logcounts(sccel.norm)[high.var.genes.scran,], logcounts(scdrop.norm)[high.var.genes.scran,])
## toc()
## 
## # extract the corrected values from MNN and combine
## merged.expr.MNN.high.var = t(Reduce(cbind, MNNcorrect.res.high.var$corrected)) # for PCA in mixOmics, need to transpose the data to have cells in rows
## #dim(merged.expr.MNN)
## 
## pca.MNN.select = mixOmics::pca(merged.expr.MNN.high.var, ncomp = 2)

## ----MNN-highvar---------------------------------------------------------
# color indicates protocol
plotIndiv(pca.MNN.select, pch = cell.line, ind.names =FALSE, group = batch, legend = TRUE, legend.title = 'Protocol', legend.title.pch = 'Cell line', title = 'MNN', xlim = c(-0.4, 0.4), ylim = c(-0.3, 0.3))


## ---- eval = eval_code---------------------------------------------------
## # tuning for MINT to decide on the number of variables to select.
## 
## # we choose 2 components here but if the interest is to discriminate all sample groups then we advise to include more components
## # see link above.
## 
## # outcome: cell line types
## Y = as.factor(cell.line[rownames(data.combined)])
## # the vector indicating each independent study / protocol or batch
## study = batch
## 
## test.keepX.values = c(seq(5,10, 5))
## 
## tune.mint = tune(X = data.combined, Y = Y, study = study, ncomp = 2, test.keepX = test.keepX.values,
## method = 'mint.splsda', dist = "max.dist", progressBar = FALSE)
## 
## # tune.mint   # lists the different types of outputs
## 
## # mean error rate per component and per tested keepX value
## # tune.mint$error.rate
## # tune.mint$choice.keepX  # the number of variables to select

## ---- eval = eval_code---------------------------------------------------
## # argument needed: how many genes to select per component:
## list.keepX = c(10,10)    # this can be an arbitrary choice, or according to our tuning parameter see code above
## # we input: tune.mint$choice.keepX
## 
## tic('MINT selection')
## mint.select = mint.splsda(X = data.combined, Y = Y, study = study, ncomp = 2, keepX = list.keepX)
## toc()

## ----MINT-select---------------------------------------------------------
# plot from mixOmics but small bug to fix with pch and style = ggplot2
#plotIndiv(mint.select, group = batch, pch = as.numeric(cell.line), pch.levels = levels(cell.line), subtitle = 'MINT with variable selection', legend =FALSE)

data.plot = data.frame(comp1 = mint.select$variates$X[,1], comp2 = mint.select$variates$X[,2])
data.plot$batch = batch
data.plot$title = "MINT"
data.plot$cell.line = cell.line


ggplot(data.plot, aes(x = comp1, y = comp2, shape = factor(cell.line))) + geom_point(aes(colour = batch), size = 3) +  scale_shape_manual(values= 1:8, name="Cell lines", labels= levels(cell.line)) +  scale_color_manual('Protocol', values = c(color.mixo(1:3)))  + labs(x = "Component 1", y = "Component 2") +  facet_grid(. ~ title) + theme_bw(base_size = 14)

## ---- eval = eval_code---------------------------------------------------
## # set up data
## sce.all.high.var = SingleCellExperiment(assays=list(counts=cbind(counts(sc10x.norm)[high.var.genes.scran,], counts(sccel.norm)[high.var.genes.scran,], counts(scdrop.norm)[high.var.genes.scran,])))
## sce.all.high.var$protocol = batch  # to accommodate for protocol effect
## sce.all.high.var$cell = cell.line  # for DE analysis
## 
## 
## tic('ZINB-WaVE high variables')
## zinb <- zinbFit(sce.all.high.var,
##                        X="~protocol",
##                        K=2,
##                        BPPARAM=BiocParallel::SerialParam())
## 
## 
## zinb.res <- zinbwave(sce.all.high.var,
##                        X="~protocol",
##                        K=2,
##                        normalizedValues=TRUE,
##                        residuals = TRUE,
##                        fitted_model = zinb,
##                        BPPARAM=BiocParallel::SerialParam(), epsilon=1e13)
## toc()
## 
## # then plot in lower dim space
## #extract latent components
## zinb.res.comp = reducedDim(zinb.res,"zinbwave")
## 
## # extract data matrix with batch removed
## weights <- assay(zinb.res, "weights")
## 

## ----ZINBWAVE-HVG--------------------------------------------------------
# plot(zinb.res.comp[, 'W1'], zinb.res.comp[, 'W2'], 
# pch = as.numeric(cell.line), col = color.mixo(as.numeric(batch)), lwd = 1.5, xlab = 'Component 1', ylab = 'Component 2', main = 'ZINB-WaVE HVG', cex.main = 2)

# plot with ggplot
data.plot = data.frame(comp1 = zinb.res.comp[, 'W1'], comp2 = zinb.res.comp[, 'W2'])
data.plot$batch = batch
data.plot$title = "ZINB-WaVe"
data.plot$cell.line = cell.line

ggplot(data.plot, aes(x = comp1, y = comp2, shape = factor(cell.line))) + geom_point(aes(colour = batch), size = 3) +  scale_shape_manual(values= 1:8, name="Cell line", labels= levels(cell.line)) +  scale_color_manual('Protocol', values=c(color.mixo(1:3)))  + labs(x = "Component 1", y = "Component 2") +  facet_grid(. ~ title) + theme_bw(base_size = 14) 

## ---- eval = eval_code---------------------------------------------------
## ncomp = 15
## # MultiCCA for more than 2 data sets
## srt_all = RunMultiCCA(list(srt10x, srtCEL, srtDROP), genes.use = high.var.genes.scran, num.ccs = ncomp)
## srt_all <- AlignSubspace(srt_all, reduction.type = "cca", grouping.var="protocol", verbose = FALSE, dims.align = 1:ncomp, num.possible.genes = length(high.var.genes.scran))

## ----Seurat-HVG----------------------------------------------------------
#t SNE plot
srt_all <- RunTSNE(object = srt_all, reduction.use = "cca.aligned", dims.use = 1:ncomp, 
    do.fast = TRUE)
p1 <- TSNEPlot(object = srt_all, group.by = "protocol", do.return = TRUE, pt.size = 0.5)
p2 <- TSNEPlot(object = srt_all, group.by = "cell line", do.return = TRUE, pt.size = 0.5)
# plot_grid(p1, p2)

# color indicates protocol
# plot(p1$data[,1], p1$data[,2], 
# pch = as.numeric(cell.line), col = color.mixo(as.numeric(batch)), lwd = 1.5, xlab = 'tSNE 1', ylab = 't-SNE 2', main = 'Seurat, high variable genes', cex.main = 2)

# with ggplot
data.plot = data.frame(comp1 = p1$data[,1], comp2 = p1$data[,2])
data.plot$batch = batch
data.plot$title = "dCCA"
data.plot$cell.line = cell.line

ggplot(data.plot, aes(x = comp1, y = comp2, shape = factor(cell.line))) + geom_point(aes(colour = batch), size = 3) +  scale_shape_manual(values= 1:8, name="Cell line", labels= levels(cell.line)) +  scale_color_manual('Protocol', values=c(color.mixo(1:3)))  + labs(x = "tSNE 1", y = "tSNE 2") +  facet_grid(. ~ title) + theme_bw(base_size = 14) 



## ---- eval = eval_code---------------------------------------------------
## sce.all = SingleCellExperiment(
##   assays=list(
##     counts=cbind(counts(sc10x.norm)[list.intersect,], counts(sccel.norm)[list.intersect,],counts(scdrop.norm)[list.intersect,]),
##     logcounts=cbind(logcounts(sc10x.norm)[list.intersect,], logcounts(sccel.norm)[list.intersect,], logcounts(scdrop.norm)[list.intersect,]))
## )
## 
## sce.all$batch = batch  # to accommodate for protocol effect
## sce.all$cell = cell.line  # for the supervised analysis, if needed
## 
## gene.var.10X = apply(assay(sc10x.norm), 1, var)
## #hist(gene.var.10X)
## 
## gene.var.drop = apply(assay(scdrop.norm), 1, var)
## #hist(gene.var.drop)
## 
## gene.var.cel = apply(assay(sccel.norm), 1, var)
## #hist(gene.var.cel)
## 
## # choose most lowly variable genes and intersection
## k = 2000
## SEG = Reduce(intersect, list(names(gene.var.10X)[1:k], names(gene.var.drop)[1:k], names(gene.var.drop)[1:k]))
## length(SEG)

## ---- eval = eval_code---------------------------------------------------
## scmerge.unsup.res <- scMerge(sce_combine = sce.all,
##                     ctl = SEG,
##                     # K can be the number of groups per study
##                     kmeansK = c(rep(nlevels(batch),3)),
##                     assay_name = "scMerge_unsupervised"
##                     )
## 
## 
## # PCA
## pca.scmerge = mixOmics::pca(t(as.matrix(scmerge.unsup.res@assays$data$scMerge_unsupervised)), ncomp = 2)

## ----scMerge-------------------------------------------------------------
# color indicates protocol
#plotIndiv(pca.scmerge, pch = cell.line, ind.names =FALSE, group = batch, legend = FALSE, legend.title = 'Protocol', legend.title.pch = 'Cell line', title = 'scMerge SEG')

# with ggplo2 for homogeneised outputs
data.plot = data.frame(comp1 = pca.scmerge$variates$X[,1], comp2 = pca.scmerge$variates$X[,2])
data.plot$batch = batch
data.plot$title = "scMerge"
data.plot$cell.line = cell.line

ggplot(data.plot, aes(x = comp1, y = comp2, shape = factor(cell.line))) + geom_point(aes(colour = batch), size = 3) +  scale_shape_manual(values= 1:8, name="Cell lines", labels= levels(cell.line)) +  scale_color_manual('Protocol', values = c(color.mixo(1:3)))  + labs(x = "Component 1", y = "Component 2") +  facet_grid(. ~ title) + theme_bw(base_size = 14)

## ---- eval = FALSE-------------------------------------------------------
## scmerge.sup.res <- scMerge(sce_combine = sce.all,
##                     ctl = SEG,
##                     kmeansK = c(3,3,3),
##                     assay_name = "scMerge_supervised",
##                     cell_type = sce.all$cell
##                     )
## 
## scmerge.sup.res <- scater::runPCA(scmerge.sup.res,
##                            exprs_values = "scMerge_supervised")
## 
## scater::plotPCA(scmerge.sup.res,
##                 colour_by = "batch",
##                 shape_by = "cell")

## ---- eval = eval_code---------------------------------------------------
## scanorama.corrected = scanorama$correct(list(t(logcounts(sc10x.norm)[high.var.genes.scran,]), t(logcounts(sccel.norm)[high.var.genes.scran,]), t(logcounts(scdrop.norm)[high.var.genes.scran,])), list(high.var.genes.scran, high.var.genes.scran, high.var.genes.scran), return_dense = TRUE)
## 
## data.scanorama = Reduce(rbind, scanorama.corrected[[1]])
## 
## pca.scanorama = mixOmics::pca(data.scanorama, ncomp = 2)

## ----scanorama-----------------------------------------------------------
# color indicates protocol
plotIndiv(pca.scanorama, pch = cell.line, ind.names =FALSE, group = batch, legend = FALSE, legend.title = 'Protocol', legend.title.pch = 'Cell line', title = 'Scanorama', xlim = c(-0.4, 0.4), ylim = c(-0.3, 0.3))

## ---- eval = eval_code---------------------------------------------------
## # on original data
## kBET.estim.orig <- kBET(t(as.matrix(data.combined)), batch, plot=FALSE, testSize = 25)
## 
## # on components
## #kBET.estim.orig.pc <- kBET(t(as.matrix(pca.combined$variates$X)), batch, plot=FALSE, testSize = 25)
## 
## #scMerge
## kBET.estim.scmerge <- kBET(t(as.matrix(scmerge.unsup.res@assays$data$scMerge_unsupervised)), batch, plot=FALSE, testSize = 25)
## 
## #kBET.estim.scmerge.pc <- kBET(t(as.matrix(pca.scmerge$variates$X)), batch, plot=FALSE, testSize = 25)
## 
## # MNN
## kBET.estim.mnn <- kBET(as.matrix(merged.expr.MNN.high.var), batch, plot=FALSE, testSize = 25)
## 
## #kBET.estim.mnn.pc <- kBET(t(as.matrix(pca.MNN.select$variates$X)), batch, plot=FALSE, testSize = 25)
## 
## #scanorama
## kBET.estim.scano <- kBET(t(as.matrix(data.scanorama)), batch, plot=FALSE, testSize = 25)
## 
## #kBET.estim.scano.pc <- kBET(t(as.matrix(pca.scanorama$variates$X)), batch, plot=FALSE, testSize = 25)
## 
## # ZINB-WaVe
## # on weights extracted
## kBET.estim.zinb <- kBET(t(as.matrix(weights)), batch, plot=TRUE, testSize = 25)
## 
## # on reduced dim
## #kBET.estim.zinb.pc <- kBET(as.matrix(zinb.res.comp), batch, plot=TRUE, testSize = 50)
## 
## # MINT
## #kBET.estim.mint.pc <- kBET(as.matrix(mint.select$variates$X), batch, plot=FALSE, testSize = 25)

## ----kBET-results--------------------------------------------------------
# data.kBET.pc = data.frame(orig.expected = kBET.estim.orig.pc$stats$kBET.expected, 
#                           orig = kBET.estim.orig.pc$stats$kBET.observed, 
#                           MNN = kBET.estim.mnn.pc$stats$kBET.observed, 
#                           ZINBWaVE = kBET.estim.zinb.pc$stats$kBET.observed, 
#                           Scanorama  = kBET.estim.scano.pc$stats$kBET.observed,
#                           scMerge = kBET.estim.scmerge.pc$stats$kBET.observed,
#                           MINT = kBET.estim.mint.pc$stats$kBET.observed)
# method = rep(colnames(data.kBET.pc), each = 100)
# method = factor(method, as.character(colnames(data.kBET.pc)))  # to reorder
# data.kBET.plot = data.frame(observed_kBET = as.vector(unlist(data.kBET.pc)), method)
# 
# # grouped boxplot
# p = ggplot(data.kBET.plot, aes(x=method, y=observed_kBET, fill = method)) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + geom_boxplot() 
# p



# on batch removed data:
data.kBET = data.frame(orig.expected = kBET.estim.orig$stats$kBET.expected, 
                      orig = kBET.estim.orig$stats$kBET.observed, 
                       MNN = kBET.estim.mnn$stats$kBET.observed, 
                       ZINBWaVE = kBET.estim.zinb$stats$kBET.observed, 
                       Scanorama  = kBET.estim.scano$stats$kBET.observed,
                       scMerge = kBET.estim.scmerge$stats$kBET.observed)
method = rep(colnames(data.kBET), each = 100)
method = factor(method, as.character(colnames(data.kBET)))  # to reorder
data.kBET.plot = data.frame(observed = as.vector(unlist(data.kBET)), method)

# grouped boxplot
p = ggplot(data.kBET.plot, aes(x=method, y=observed, fill = method)) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + geom_hline(yintercept=0.75, linetype="dashed", color = "lightblue") +  geom_boxplot() 
p


## ------------------------------------------------------------------------
# function that calculates the silhouette coefficient based on a known cluster (i.e. protocol - batch or cell line)
# calculates silhouette width average
calc.sil = function(
  x, # the PC variates
  y1, y2 = NULL, # factor of interest, e.g. known batch info or known cell type
  name.y1, name.y2 = NULL # character of the factor of interest
){
  library(cluster)
  # calculate the distance, here euclidean is appropriate for PCA, NOT for t-SNE
  dist.res = daisy(x, metric = 'euclidean')
  # for factor 1
  sil.batch.res1 = silhouette(x = as.numeric(y1), dist = dist.res)
  # if factor 2 is provided
  if(!is.null(y2))  sil.batch.res2 = silhouette(x = as.numeric(y2), dist = dist.res)
  
  # extract average width silhouette per level
  res1 = c(summary(sil.batch.res1)["clus.avg.widths"]$clus.avg.widths)
  names(res1) = levels(y1)
  if(!is.null(y2)){
    res2 = c(summary(sil.batch.res2)["clus.avg.widths"]$clus.avg.widths)
    names(res2) = levels(y2)
  }
  
  # output data for plotting
  if(!is.null(y2)){
    silh.coeff = c(res1, res2)
    cluster = c(levels(y1), levels (y2))
    type = c(rep(name.y1, nlevels(y1)), rep(name.y2, nlevels(y2)))
  }else{
    silh.coeff = c(res1)
    cluster = c(levels(y1))
    type = rep(name.y1, nlevels(y1))
  }

  data.plot = data.frame(silh.coeff, cluster, type)

  
  return(invisible(data.plot))
}


## ----silhouette----------------------------------------------------------
## original data
silh.orig = calc.sil(x = pca.combined$variates$X, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

# MNN
silh.MNN = calc.sil(x = pca.MNN.select$variates$X, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

# MINT
silh.MINT = calc.sil(x = mint.select$variates$X,y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

# zinb-wave
silh.ZINB = calc.sil(x = zinb.res.comp, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

# scMerge
silh.scMerge = calc.sil(x = pca.scmerge$variates$X, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

# scanorama
silh.scano = calc.sil(x = pca.scanorama$variates$X, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')


# merge all results for plotting
data.plot = rbind(silh.orig, silh.MNN, silh.MINT, silh.ZINB, silh.scMerge, silh.scano)
data.plot$method = c(rep('orig', nrow(silh.orig)), 
                     rep('MNN', nrow(silh.MNN)),
                     rep('MINT', nrow(silh.MINT)),
                     rep('ZINB', nrow(silh.ZINB)),
                     rep('scMerge', nrow(silh.scMerge)),
                     rep('scanorama', nrow(silh.scano))
)
data.plot$method = factor(data.plot$method, levels = unique(data.plot$method)) # to reorder

ggplot(data.plot, aes(x=type, y=silh.coeff, fill = type)) + geom_boxplot() + facet_grid(cols = vars(method)) + theme(axis.text.x = element_text(angle = 60, hjust = 1), strip.text = element_text(size=10)) + scale_fill_manual(values=c("#999999", "#E69F00")) + labs(x = "Cluster type", y = "Silhouette Coefficient", name="Cluster type") 


## ------------------------------------------------------------------------
calc.ARI = function(
  x, # the PC variates
  y1, y2 = NULL, # factor of interest, e.g. known batch info or known cell type
  name.y1, name.y2 = NULL # character of the factor of interest
){
  library(clues)
  library(cluster)
  # calculate the distance, here euclidean is appropriate for PCA, NOT for t-SNE
  dist.res = daisy(x, metric = 'euclidean')
  
  # need to cluster the data, eg. from pam
  # for factor 1
  pam.res1 = pam(dist.res, diss = TRUE, k = nlevels(y1))
  res1 = adjustedRand(pam.res1$clustering, as.numeric(y1))
  # for factor 2 if provided
  if(!is.null(y2)){
    pam.res2 = pam(dist.res, diss = TRUE, k = nlevels(y2))
    res2 = adjustedRand(pam.res2$clustering, as.numeric(y2))
  }
  
  if(!is.null(y2)){
     res = rbind(res1, res2) 
     rownames(res) = c(name.y1, name.y2)
  }else{
    res = res1
    rownames(res) = name.y1
    }
  
  return(invisible(res))
}

## ------------------------------------------------------------------------
## original data
ari.orig = calc.ARI(x = pca.combined$variates$X, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

# MNN
ari.MNN = calc.ARI(x = pca.MNN.select$variates$X, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

# MINT
ari.MINT = calc.ARI(x = mint.select$variates$X, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

# zinb-wave
ari.ZINB = calc.ARI(x = zinb.res.comp, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

# scMerge
ari.scMerge = calc.ARI(x = pca.scmerge$variates$X, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

# scanorama
ari.scano = calc.ARI(x = pca.scanorama$variates$X, y1 = batch, y2 = cell.line, name.y1 = 'batch', name.y2 = 'cell line')

#output
data.ARI = rbind(ari.orig, ari.MNN, ari.MINT, ari.ZINB, ari.scMerge, ari.scano)
method = c(rep('orig', nrow(ari.orig)), 
                     rep('MNN', nrow(ari.MNN)),
                     rep('MINT', nrow(ari.MINT)),
                     rep('ZINB', nrow(ari.ZINB)),
                     rep('scMerge', nrow(ari.scMerge)),
                     rep('scanorama', nrow(ari.scano))
)


# kable
dt = rbind(method, ARI = round(data.ARI[, 'Rand'],2))
kable(t(dt), caption = 'Summary ARI for components-based methods', digits = 2) %>%
  kable_styling(bootstrap_options = "striped", font_size = 7)

write.csv(t(dt), 'ARI-puremix.res.csv')

## ------------------------------------------------------------------------
sessionInfo()

## ---- echo = FALSE, eval = FALSE-----------------------------------------
## # render as R document
## purl('CellBench_puremix_integration.Rmd')

## ---- echo = FALSE, eval = eval_code-------------------------------------
## save.image('CellBench_puremix.results.RData')

