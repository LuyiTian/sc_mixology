# Cellbench

[CellBench](https://github.com/LuyiTian/CellBench_data) uses three human lung adenocarcinoma cell lines HCC827, H1975 and H2228, which were cultured separately, and the same batch was processed in three different ways. Firstly, the single cells were mixed in equal proportions, with libraries generated using three different protocols: CEL-seq2, Drop-seq with Dolomite equipment and 10X Chromium. Secondly, the single cells were sorted from the three cell lines into 384-well plates, with an equal number of cells per well in different combinations. Thirdly, the RNA were extracted in bulk for each cell line and the RNA was mixed in 7 different proportions and diluted to single cell equivalent amounts.

<img src=script/expr_design.png width="800">

The dataset is stored in R object using SingleCellExperiment class. Below are instructions for getting four files: metadata (including annotations) and count data for each dataset. All data was after quality control, without gene filtering.


## read files for R

You can find R object files in the [data](https://github.com/LuyiTian/CellBench_data/tree/master/data) folder 

```R
load("data/sincell_with_class.RData")
```

## Counts

As specified by SingleCellExperiment, the count can be retrived by  `counts(sce)` function:

```R
counts(sce10x_qc)[1:5, 1:5]
```

## Metadata

As specified by SingleCellExperiment, the Metadata can be retrived by  `colData(sce)` function:

```R
head(colData(sce10x_qc))
```


## examples of using these dataset

You can find the Rnotebook on [script](https://github.com/LuyiTian/CellBench_data/tree/master/script) folder. `data_explore_mixture.Rmd` includes codes for analysing cell mixture and RNA mixture dataset. `data_integration_single_cell.Rmd` shows example for data integration for single cell dataset from three different protocols.







The original data will be on GEO soon.

