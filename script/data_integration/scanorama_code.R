# run scanorama
# scanorama cannot be installed on linux server, so it was ran separately on a mac laptop.
library(CellBench)
library(reticulate)
library(scater)
library(scran)
log_file =paste("log_file/data_int_scanorama_all",format(Sys.time(), "%a_%b_%d"),"txt",sep = ".")
cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "start preprocessing...\n"), file = log_file, append = TRUE)


scanorama <- import('scanorama')
setwd("~/Dropbox/research/benchmark/resubmit_code")

# sce = RNAmix_data_int_raw$result[[24]]
# sce = sce[rowData(sce)$hi_var,]
# sce_list = lapply(unique(sce$batch),function(ba){t(logcounts(sce[,sce$batch==ba]))})
# 
# ge_list = lapply(unique(sce$batch),function(x){rownames(sce)[rowData(sce)$hi_var]}) # they are the same
# 
# scanorama.corrected = scanorama$correct(sce_list, ge_list, return_dense = TRUE)
# 
# data.scanorama = Reduce(rbind, scanorama.corrected[[1]])
# 
# 
# 
# assay(sce,"int_expr") = t(data.scanorama)
# 
# plotPCA(sce,run_args=list(exprs_values = "logcounts"),colour_by="batch")
# 
# plotPCA(sce,run_args=list(exprs_values = "int_expr"),colour_by="batch")



scanorama_int = function(sce){
  tp = system.time({
    try_res = try({
      if("logcounts" %in% assayNames(sce)){
        if (max(logcounts(sce))>100){
          logcounts(sce) = log2(logcounts(sce)+1)
        }
        sce = sce[rowData(sce)$hi_var,]
        sce_list = lapply(unique(sce$batch),function(ba){t(logcounts(sce[,sce$batch==ba]))})
        ge_list = lapply(unique(sce$batch),function(x){rownames(sce)[rowData(sce)$hi_var]}) # they are the same
        scanorama.corrected = scanorama$correct(sce_list, ge_list, return_dense = TRUE)
        data.scanorama = Reduce(rbind, scanorama.corrected[[1]])
        assay(sce,"int_expr") = t(data.scanorama)        
      }
    })
    if (class(try_res) == "try-error") {
      cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
    }
  })
  
  method_name = "scanorama"
  method_type = "data_integration"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}



data_int_method <- list(
  scanorama=scanorama_int
)



cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "start to apply data integration methods...\n"), file = log_file, append = TRUE)

for (ix in 9){
  print(ix)
  res = NULL
  res <- readRDS(paste0("~/data/rdata/RNAmix_int_data_",ix,".Rds"))
  res_t <- res %>%
    apply_methods(data_int_method)
  saveRDS(res_t, file=paste0("~/data/rdata/RNAmix_after_scanorama_",ix,".Rds"))
  res_t = NULL
}


for (ix in 1:9){
  print(ix)
  res = NULL
  res <- readRDS(paste0("~/data/rdata/sc_int_data_",ix,".Rds"))
  res_t <- res %>%
    apply_methods(data_int_method)
  saveRDS(res_t, file=paste0("~/data/rdata/sc_after_scanorama_",ix,".Rds"))
  res_t = NULL
}



cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "finish all data integration methods...\n"), file = log_file, append = TRUE)





