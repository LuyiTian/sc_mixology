## data integration
library(scater)
library(scran)
library(CellBench)
library(scPipe)
library(R.utils)
MAX_TIME = 60*120
setwd("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit")
log_file =paste("log_file/data_int_SC",format(Sys.time(), "%a_%b_%d"),"txt",sep = ".")
cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "start preprocessing...\n"), file = log_file, append = TRUE)

res2 = readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/final/sc_all_after_imputation.Rds")

# feature selection
scran_high_var = function(sce,topn=1500){ # to have more than 1000 so after combine the hi_var genes will be close to 1000
  var.fit <- trendVar(sce, method="loess", use.spikes=FALSE)
  var.out <- decomposeVar(sce, var.fit)
  hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
  rowData(sce)$hi_var = FALSE
  rowData(sce)$hi_var[rownames(rowData(sce)) %in% rownames(hvg.out)] = TRUE
  return(sce)
}

scran_low_var = function(sce,topn=1000){ # for scMerge
  var.fit <- trendVar(sce, method="loess", use.spikes=FALSE)
  var.out <- decomposeVar(sce, var.fit)
  var.out = var.out[var.out$mean>mean(var.out$mean),]
  lvg.out <- var.out[order(var.out$bio, decreasing=FALSE)[1:topn], ]
  rowData(sce)$lo_var = FALSE
  rowData(sce)$lo_var[rownames(rowData(sce)) %in% rownames(lvg.out)] = TRUE
  return(sce)
}

hivar_method = list(scran_hi = scran_high_var)
res2 = res2 %>% apply_methods(hivar_method)
lovar_method = list(scran_hi = scran_low_var)
res2 = res2 %>% apply_methods(lovar_method)

res2 = res2[,!(colnames(res2) %in% c("hivar_method","lovar_method"))]
res2 = res2[!(res2$data %in% c("sc_Celseq2_5cl_p1" ,"sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")),]
res2 = res2[!(res2$norm_method =="BASiCS"),] # BASiCS does not work on droplet based protocol.


cmb_sce = function(sce_list,batch_names){
  for (ix in 1:length(sce_list)){ # to have consist gene id type.
    if(gene_id_type(sce_list[[ix]])=="external_gene_name"){
      rownames(sce_list[[ix]]) = sce_list[[ix]]@int_elementMetadata[, "ensembl_gene_id"]
      gene_id_type(sce_list[[ix]])="ensembl_gene_id"
    }
  }
  com_genes = Reduce(intersect,lapply(sce_list,rownames))
  com_colna = Reduce(intersect,lapply(sce_list,function(x){colnames(colData(x))}))
  com_hi_var = Reduce(and,lapply(sce_list,function(x){rowData(x[com_genes,])$hi_var}))
  com_lo_var = Reduce(and,lapply(sce_list,function(x){rowData(x[com_genes,])$lo_var}))
  for (ix in 1:length(sce_list)){
    sce_list[[ix]] = sce_list[[ix]][com_genes,]
    colData(sce_list[[ix]]) = colData(sce_list[[ix]])[,com_colna]
    colnames(sce_list[[ix]]) = paste(batch_names[ix], colnames(sce_list[[ix]]), sep="_")
    assays(sce_list[[ix]]) = list(counts=counts(sce_list[[ix]]),logcounts=logcounts(sce_list[[ix]]))
    if(!("batch" %in% colnames(colData(sce_list[[ix]])))){
      sce_list[[ix]]$batch = batch_names[ix]
      rowData(sce_list[[ix]])$hi_var = com_hi_var
      rowData(sce_list[[ix]])$lo_var = com_lo_var # for scMerge
    }
  }
  combined_sce = Reduce(cbind, sce_list)
  return(combined_sce)
}

res_int = res2[1:2,]

for (vnorm in as.character(unique(res2$norm_method))) {
  for(vimp in as.character(unique(res2$impute_method))){
    print(vnorm)
    print(vimp)
    tmp = res2[res2$norm_method==vnorm & res2$impute_method==vimp,]
    if(all(unlist(lapply(tmp$result,function(x){all(c("counts","logcounts") %in% assayNames(x))})))){
      c_sce = cmb_sce(tmp$result, as.character(tmp$data))
      res_df = data.frame(data="single_cell",norm_method=vnorm,impute_method=vimp)
      res_df$result <- list(c_sce)
      res_int = rbind(res_int,res_df)
    }
  }
}
res_int = res_int[-(1:2),]

for(i in 1:8){
  saveRDS(res_int[(4*(i-1)+1):(4*i),],paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/sc_int_data_",i ,".Rds"))
}



################################
################################
library(Seurat)
library(scMerge)
library(zinbwave)
library(mixOmics)

no_int = function(sce){
  tp = system.time({
    if (max(logcounts(sce))>100){
      logcounts(sce) = log2(logcounts(sce)+1)
    }
    assay(sce,"int_expr") = logcounts(sce)
  })
  method_name = "no_integration"
  method_type = "data_integration"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}


seurat_int = function(sce){
  tp = system.time({
    try_res = try({
      if (max(logcounts(sce))>100){
        logcounts(sce) = log2(logcounts(sce)+1)
      }
      srt = CreateSeuratObject(raw.data = counts(sce),display.progress = FALSE)
      srt@scale.data = logcounts(sce)
      if("hi_var" %in% colnames(rowData(sce))){
        srt@var.genes = rownames(sce)[rowData(sce)$hi_var]
      }
      srt@meta.data[, "batch"] <- sce$batch
      srt_list = lapply(unique(sce$batch),function(ba){SubsetData(srt,cells.use = srt@cell.names[srt@meta.data[, "batch"]==ba])})
      if(length(srt_list)==2){
        srt_comb = RunCCA(srt_list[[1]],srt_list[[2]],genes.use=srt@var.genes,num.cc=10)
      }else{
        srt_comb = RunMultiCCA(srt_list,genes.use=srt@var.genes,num.ccs=10)
      }
      srt_comb = AlignSubspace(srt_comb, reduction.type = "cca", grouping.var="batch", verbose = FALSE, dims.align = 1:10, num.possible.genes = length(srt@var.genes))
      srt_comb = RunTSNE(object = srt_comb, reduction.use = "cca.aligned", dims.use = 1:10, 
                         do.fast = TRUE,check_duplicates = FALSE)
      reducedDim(sce,"TSNE") = srt_comb@dr$tsne@cell.embeddings
    })
    if (class(try_res) == "try-error") {
      cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
    }
  })
  
  method_name = "Seurat"
  method_type = "data_integration"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}



mnn_int = function(sce){
  tp = system.time({
    try_res = try({
      withTimeout({
        if (max(logcounts(sce))>100){
          logcounts(sce) = log2(logcounts(sce)+1)
        }
        sce = sce[rowData(sce)$hi_var,]
        lc_list = lapply(unique(sce$batch),function(ba){logcounts(sce[,sce$batch==ba])})
        sce_corrected = do.call(mnnCorrect,lc_list)
        comb_mtx = Reduce(cbind,sce_corrected$corrected)
        colnames(comb_mtx) = colnames(sce)
        assay(sce,"int_expr")  = comb_mtx
        
      },timeout = MAX_TIME, onTimeout="error")
    })
    if (class(try_res) == "try-error") {
      cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
    }
  })
  
  method_name = "MNNs"
  method_type = "data_integration"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}


scMerge_int = function(sce){
  tp = system.time({
    try_res = try({
      if (max(logcounts(sce))>100){
        logcounts(sce) = log2(logcounts(sce)+1)
      }
      sce = scMerge(sce_combine = sce, 
                    ctl =rownames(sce)[rowData(sce)$lo_var],
                    marker=  rownames(sce)[rowData(sce)$hi_var],
                    kmeansK = rep(10,length(unique(sce$batch))),
                    assay_name = "int_expr")
    })
    if (class(try_res) == "try-error") {
      cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
    }
  })
  
  method_name = "scMerge_us"
  method_type = "data_integration"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}


scMerge_s_int = function(sce){
  tp = system.time({
    try_res = try({
      if (max(logcounts(sce))>100){
        logcounts(sce) = log2(logcounts(sce)+1)
      }
      if("group" %in% colnames(colData(sce))){
        ct = sce$group
      }else{
        ct = sce$cell_line_demuxlet
      }
      sce = scMerge(sce_combine = sce, 
                    ctl =rownames(sce)[rowData(sce)$lo_var],
                    marker=  rownames(sce)[rowData(sce)$hi_var],
                    kmeansK = rep(10,length(unique(sce$batch))),
                    cell_type=ct,
                    assay_name = "int_expr")
    })
    if (class(try_res) == "try-error") {
      cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
    }
  })
  
  method_name = "scMerge_s"
  method_type = "data_integration"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}


zinbwave_int = function(sce){
  tp = system.time({
    try_res = try({
      if (max(logcounts(sce))>100){
        logcounts(sce) = log2(logcounts(sce)+1)
      }
      zinb.res <- zinbwave(sce[rowData(sce)$hi_var,], 
                           X="~batch",
                           K=2,
                           normalizedValues=FALSE,
                           residuals = FALSE, epsilon=1e13,
                           BPPARAM=MulticoreParam(NUM_OF_THREAD))
      sce = sce[rownames(zinb.res),colnames(zinb.res)]
      reducedDim(sce, "zinbwave") = reducedDim(zinb.res, "zinbwave")
    })
    if (class(try_res) == "try-error") {
      cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
    }
  })
  
  method_name = "zinbwave"
  method_type = "data_integration"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}


mixOmics_mint = function(sce = sce, ## a combined sce with batch info in sce$batch 
                         #colData.class = "mix", ## name of the colData that includes the biological groups
                         use.hvgs = F, ## whether to use HVGs in rowData(sce)$hi_vars (T) or all genes (F)
                         ######### below can be left to default
                         ncomp = 5L, ## integer: good rule of thumb: number of biological groups - 1
                         keepX = c(50,40,30,20,10), ## number of genes to select on each of the ncomp components
                         tune.hps = F, ## whether to tune the number of genes on each component (T) or use keepX (F)
                         tune.keepX = seq(10,100,10), ##  if tune.hps=T: span of number of genes to explore for tuning
                         optimum.ncomp = F, ## whether to find the optimum number of sPLSDA components (T) or just use ncomp (F)
                         dist = "max.dist", ## if tune.hps: distance(s) to use for classification error rate.
                         ## subset of c("min.distance", "centroids.dist","mahalanobis.dist" ); or "all"
                         
                         measure = "BER" # if tune.hps=T: measure to use for classification error rate; one of c("overall","BER")
){
  tp = system.time({
    try_res = try({
      
      ######################### the t.test.process fuction for optimisation
      t.test.process = function(mat.error.rate, alpha = 0.01)
      {
        # mat.error.rate has nrep rows and ncomp columns
        # we test successively whether adding a component improves the results
        
        max = ncol(mat.error.rate) #number max of components included
        pval = NULL
        opt = 1 #initialise the first optimal number of components
        for(opt in 1:max)
        {
          j=opt+1
          temp = try(t.test(mat.error.rate[,opt],mat.error.rate[,j],alternative="greater")$p.value, silent=T) #t.test of "is adding X comp improves the overall results"
          if(any(class(temp) == "try-error") || is.na(temp)) # temp can be NaN when error.keepX is constant
          {
            pval = 1
          } else {
            pval = temp
          }
          #print(opt)
          #print(j)
          #print(pval)
          
          while(pval> (alpha) & j<max)
          {
            j=j+1
            temp = try(t.test(mat.error.rate[,opt],mat.error.rate[,j],alternative="greater")$p.value, silent=T) #t.test of "is adding X comp improves the overall results"
            if(any(class(temp) == "try-error") || is.na(temp)) # temp can be NaN when error.keepX is constant
            {
              pval = 1
            } else {
              pval = temp
            }
          }
          
          if( (pval> (alpha))) #if all pvalues were greater than alpha, then we do not increase opt and get out of the loop
            break
        }
        ncomp_opt = opt
        
        return(ncomp_opt)
      }
      ######################### < entry checks >
      if("group" %in% colnames(colData(sce))){
        Y = as.factor(sce$group) ## cell types/classes
      }else{
        Y = as.factor(sce$cell_line_demuxlet) ## cell types/classes
      }
      
      batch = as.factor(sce$batch) ## batch vector
      {
        ## ncomp and keepX match - otherwise all genes will be included in excess comps
        if (length(keepX)<ncomp)
          stop ("The length of keepX should be ncomp")
        
        ## if ncomp optimisation is required, ncomp>1
        if(optimum.ncomp & !isTRUE(ncomp>1))
          stop("ncomp must be an integer greater than 1 for optimisation")
        
        ## sce checks
        if(class(sce)!="SingleCellExperiment") ## make sure it is SCE object
          stop("sce must be a SingleCellExperiment object with batch metadata")
        
        ## batch checks
        if(is.null(batch)){
          stop("$batch data do not exist in rowData(sce)")
        } else if(!isTRUE(length(unique(batch))>1)){
          stop("there must be more than one batch in the data to use mint.splsda")
        }
        
        ## cell class checks
        if(is.null(Y)){
          stop("cell type data do not exist in colData(sce)")
        } else if(!isTRUE(length(unique(Y))>1)){
          stop("there must be more than one cell type in the data to use mint.splsda")
        }
        
        ## logicals
        if(!all(is.logical(use.hvgs),is.logical(tune.hps), is.logical(optimum.ncomp)))
          stop("use.hvgs, tune.hps and optimum.ncomp must be logical")
        
        ## ensure there are no duplicate cell names
        if(any(colnames(sce)!= make.unique(colnames(sce)))){
          message("changed duplicate cell names")
          colnames(sce) = make.unique(colnames(sce))
        }
        
        ## hi_var exist
        if(use.hvgs&is.null(rowData(sce)$hi_var)){
          message("HVGs are not specified in rowData(sce)$hi_var - continuing with all genes")
        }
        
        if(optimum.ncomp & nlevels(batch)<3)
          warning("The number of components cannot be reliably optimised
                  since there are less than 3 batches. Regard the results with care.
                  Refer to mixOmics documentation for details.")
        
        ## MINT checker checks the rest
        
        ######################### < /entry checks >
        
        if (max(logcounts(sce))>100){
          logcounts(sce) = log2(logcounts(sce)+1)
        }
        
        ## reduce sce if only HVGs are needed
        if(use.hvgs){
          sce = sce[rowData(sce)$hi_var,]
        }
      }
      ######################################
      ################## < wrapper > #######
      ######################################
      
      ## check if it need to be tuned or optimised:
      
      if(tune.hps){
        ## tune the number of markers
        mint.tune = tune.mint.splsda(
          X = t(logcounts(sce)),
          Y = Y,
          study =  batch,
          ## get the optimum for measure and dist
          ncomp=ncomp,
          ## assess tune.keepX number of features:
          test.keepX = tune.keepX,
          ## use dist to estimate the classification error rate
          dist = dist,
          measure=measure,
          progressBar = F)
        
        ## change keepX to the tuned vector
        keepX = mint.tune$choice.keepX
        ## choose the optimum ncomp if required
        if(optimum.ncomp){
          ncomp = t.test.process(mint.tune$error.rate)
        }
      }
      
      if(optimum.ncomp&!tune.hps){ ## if only ncomp is to be optimised
        
        ## get the tuned unoptimised sparse MINT
        mint.res= mint.splsda(X=t(logcounts(sce)),
                              Y = Y,
                              study = batch,
                              ncomp = ncomp,
                              keepX = keepX
        )
        
        ## evaluate method's performance at different no. of comp.s
        mint.performance = perf(mint.res, progressBar = F)
        ncomp = t.test.process(mint.performance$global.error$error.rate.class[[dist]])
      }
      
      
      ## final model using tuned parameters
      mint.res = mint.splsda(
        X = t(logcounts(sce)),
        Y = Y,
        study = batch,
        ncomp=ncomp,
        keepX = keepX
      )
      
      ## marker genes
      comp =1
      markers = NULL
      while(comp <= ncomp){
        markers = c(markers, selectVar(mint.res, comp = comp)$name)
        comp = comp+1
      }
      ## add a logical rowData as to whether the gene is a marker
      rowData(sce)$mint.markers <- row.names(sce) %in% markers
      ## add the sPLSDA variates for visualisation
      reducedDim(sce, "mint_variates") = mint.res$variates$X
      
      ######################################
      ################## < /wrapper > #######
      ######################################
      
    })
    if (class(try_res) == "try-error") {
      cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
    }
  })
  
  method_name = "mixOmics_mint"
  method_type = "data_integration"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}



data_int_method <- list(
  no_int=no_int,
  Seurat=seurat_int,
  scMerge_s=scMerge_s_int,
  #zinbwave=zinbwave_int,
  scMerge_us=scMerge_int,
  MNNs=mnn_int,
  mixOmics_mint=mixOmics_mint
)



cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "start to apply data integration methods...\n"), file = log_file, append = TRUE)

for (i in 1:8){
  res_int = NULL
  res_int = readRDS(paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/sc_int_data_",i ,".Rds"))
  data_int_method <- list(
    no_int=no_int,
    Seurat=seurat_int,
    scMerge_s=scMerge_s_int,
    #zinbwave=zinbwave_int,
    scMerge_us=scMerge_int,
    MNNs=mnn_int,
    mixOmics_mint=mixOmics_mint
  )
  res_t <- res_int %>%
    apply_methods(data_int_method)
  data_int_method = list(zinbwave=zinbwave_int)
  res_t1 = res_int[res_int$norm_method=="none" & res_int$impute_method=="no_impute", ] %>%
    apply_methods(data_int_method)
  res_t = rbind(res_t, res_t1)
    
    saveRDS(res_t, file=paste0("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/analysis_for_resubmit/rdata/data_int/sc_all_after_data_int_",i ,".Rds"))
    res_t = NULL
}


cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "finish all data integration methods...\n"), file = log_file, append = TRUE)



