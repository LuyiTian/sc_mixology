library(SingleCellExperiment)
library(scMerge) ## to combine sce's
library(mixOmics)

######################### < INPUT >
## load individual sce object
## if they're combined already, ensure the colData match
load("sincell_with_class.RData")
## create a list of sce's - must be clearly named by batch
sce.list = list("10X" = sce_sc_10x_qc,
                "CELseq2" = sce_sc_CELseq2_qc,
                "Dropseq" = sce_sc_Dropseq_qc)
cell.line.name = "cell_line" ## the colData for biological group of cells, cell_line, mix etc
######################### < /INPUT >



######################### combine
## combine sce's
sce.sincell =
  scMerge::sce_cbind(
    sce.list,method = "intersect", ## instersect genes
    batch_names = names(sce.list), ## include batch names
    colData_names = cell.line.name, ## include cell lines
    cut_off_batch = 0, ## ?sce_cbind
    cut_off_overall = 0) ## ?sce_cbind

######################### wrapper

mint.wrapper = function(sce = sce.sincell, ## a combined sce with sce$batch that
                        ## icludes $batch and $cell_line (or similar)
                        ncomp = 4, ## no need to change this - could put to 4 for initial runs
                        keepX = 30,
                        tune.keepX = c(5,10,35,seq(100,1000,100))
                        ){
  
  stopifnot(all(
    length(unique(sce$batch))>1, ## ensure various batches exist
    length(unique(sce[[cell.line.name]]))>1 ## ensure various cell lines exist
  ))
  ## ensure there are no duplicate cell names
  colnames(sce) = make.unique(colnames(sce))
  
  ## get the untuned unoptimised sparse MINT
  mint.prelim = mint.splsda(X=t(logcounts(sce)),
                            Y = sce$cell_line,
                            study = sce$batch,
                            ncomp = ncomp,
                            keepX = rep(keepX, ncomp)
                              )
  ## evaluate method's performance at different no. of comp.s
  mint.performance = perf(mint.prelim, progressBar = F)

  ## tune the number of markers
  tune.mint = tune.mint.splsda(
    X = t(logcounts(sce)),
    Y = sce$cell_line,
    study = factor(sce$batch),
    ## get the optimum for BER and max.distance
    ncomp=mint.performance$choice.ncomp[2,1],
    ## assess numbers 5, 35,100,500,1000:
    test.keepX = tune.keepX,
    ## use all distances to estimate the classification error rate
    dist = c('max.dist'),
    progressBar = F
  )
  
  ## final model using tuned parameters
  mint.splsda.tuned = mint.splsda(
    X = t(logcounts(sce)),
    Y = sce$cell_line,
    study = sce$batch,
    ncomp=mint.performance$choice.ncomp[2,1],
    keepX = tune.mint$choice.keepX
  )
  
  ## marker genes
  comp =1
  markers = NULL
  while(comp <= mint.performance$choice.ncomp[2,1]){
    markers = c(markers, selectVar(mint.splsda.tuned, comp = comp)$name)
    comp = comp+1
  }
  # ## plot the tuned mint.splsda plot for the combined dataset
  # plotIndiv(mint.out$mint.splsda.tuned, study = 'global', legend = T,
  #           title = 'MINT sPLS-DA',  subtitle = 'Global', ellipse=T, legend.title = 'Cell Line')
  # 
 return(list(mint.object = mint.splsda.tuned , markers = markers))
}


######################### run
mint.out = mint.wrapper(sce = sce.sincell)


## sample plot
## plot the tuned mint.splsda plot for the combined dataset
plotIndiv(mint.out$mint.object, study = 'global', legend = T,
          title = 'MINT sPLS-DA',  subtitle = 'Global', ellipse=T, legend.title = 'Cell Line')

## markers
mint.out$markers