library(BiocParallel)

setwd("~/grpu_mritchie_1/LuyiTian/SCmixology_code/benchmark_cluster/")


output_link_1 <- "Clustering_Results/"
algorithm_link <- "Clustering_Algorithms/"
input_data_folder <- "Data/"
final_output <- "Clustering_Result.RData"

input_data <- list.files(input_data_folder)[5:6]
input_data <- paste0(input_data_folder, input_data)

algorithm <- c("RaceID", "RaceID2", "RCA", "sc3", "Seurat", "clusterExperiment")

for (i in 1:length(input_data)){

  output_link <- paste0(output_link_1, gsub(".RData", "", strsplit(input_data[i], "/")[[1]][2]), "/")
  system(paste0("mkdir ", output_link))
  arguments <- list(output_link = output_link, algorithm_link = algorithm_link, algorithm = algorithm, input_data = input_data[i])

  call_algo <- function(x, arguments) {
      system.time(system( paste0( "Rscript --vanilla " , arguments$algorithm_link , arguments$algorithm[x], "_call.R ", arguments$input_data,
          " ", arguments$output_link, arguments$algorithm[x], ".txt") ) )
  }

  multicoreParam <- MulticoreParam(workers = 6)
  time_taken <- bplapply(1:6, call_algo, arguments=arguments,  BPPARAM = multicoreParam)
  names(time_taken) <- algorithm
  
  save(time_taken, file=paste0(output_link, "Time_Taken.RData"))
  
}
  



