# Cellbench

[CellBench](https://github.com/LuyiTian/CellBench_data) uses three human lung adenocarcinoma cell lines HCC827, H1975 and H2228, which were cultured separately, and the same batch was processed in three different ways. Firstly, the single cells were mixed in equal proportions, with libraries generated using three different protocols: CEL-seq2, Drop-seq with Dolomite equipment and 10X Chromium. Secondly, the single cells were sorted from the three cell lines into 384-well plates, with an equal number of cells per well in different combinations. Thirdly, the RNA was extracted in bulk for each cell line and the RNA was mixed in 7 different proportions and diluted to single cell equivalent amounts.

<img src=experiment_design/expr_design.png width="800">

The dataset is stored in R object using SingleCellExperiment class. Below are instructions for getting four files: metadata (including annotations) and count data for each dataset. All data was after quality control, without gene filtering.

## summary of all dataset

<img src=experiment_design/supp_table.png width="800">

<img src=experiment_design/supp_table_design.png width="800">


## read files for R

You can find R object files in the [data](https://github.com/LuyiTian/CellBench_data/tree/master/data) folder 

```R
load("data/sincell_with_class.RData")
```

## Counts

As specified by SingleCellExperiment, the count can be retrieved by  `counts(sce)` function:

```R
counts(sce10x_qc)[1:5, 1:5]
```

## Metadata

As specified by SingleCellExperiment, the Metadata can be retrieved by  `colData(sce)` function:

```R
head(colData(sce10x_qc))
```


## examples of using these dataset

You can find the Rnotebook on [script/data_QC_visualization](https://github.com/LuyiTian/CellBench_data/tree/master/script/data_QC_visualization) folder. `data_explore_mixture.Rmd` includes codes for analysing cell mixture and RNA mixture dataset.

## script for reproducing the analysis

the [script] folder contains scripts that can reproduce the analysis and figures in paper: scRNA-seq mixology: towards better benchmarking of single cell RNA-seq protocols and analysis methods. 

***Note:*** the `ggtern` package, which has been used to generate the ternary plot, has known issue with the recent version of `ggplot` and the relevant code is broken if you have updated `ggplot` package.



