# CellBench: single cell RNA-seq benchmarking

[CellBench](https://github.com/LuyiTian/CellBench_data) uses three human lung adenocarcinoma cell lines HCC827, H1975 and H2228, which were cultured separately, and then processed in three different ways. Firstly, single cells from each cell line were mixed in equal proportions, with libraries generated using three different protocols: CEL-seq2, Drop-seq (with Dolomite equipment) and 10X Chromium. Secondly, the single cells were sorted from the three cell lines into 384-well plates, with an equal number of cells per well in different combinations (generally 9-cells, but with some 90-cell `population` controls). Thirdly, RNA was extracted in bulk for each cell line and the RNA was mixed in 7 different proportions and diluted to single cell equivalent amounts ranging from 3.75pg to 30pg and preocessed using CEL-seq 2 and SORT-seq. ERCC spike-in controls were present in samples processed using the 2 plate-based technologies (CEL-seq2 and SORT-seq).

<img src=experiment_design/expr_design.png width="800">

Raw data from this series of experiments is available under GEO accession number [GSE118767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118767).
The processed count data obtained from [scPipe](https://bioconductor.org/packages/release/bioc/html/scPipe.html) is stored in R objects that use the [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) class. Below are instructions for getting the count data and metadata (including annotations) for each dataset. All data is post sample quality control, without gene filtering.

## Summary of all datasets

<img src=experiment_design/supp_table.png width="800">

<img src=experiment_design/supp_table_design.png width="800">

## Load files into R

You can find R object files in the [data](https://github.com/LuyiTian/CellBench_data/tree/master/data) folder 

```R
load("data/sincell_with_class.RData")
```

## Counts

To access count data from a SingleCellExperiment object, use the `counts(sce)` function:

```R
counts(sce10x_qc)[1:5, 1:5]
```

## Metadata

To access sample information from a SingleCellExperiment object, use the `colData(sce)` function:

```R
head(colData(sce10x_qc))
```

## Examples of using these datasets

You can find an Rnotebook in the [script/data_QC_visualization](https://github.com/LuyiTian/CellBench_data/tree/master/script/data_QC_visualization) folder named `data_explore_mixture.Rmd` which includes code for analysing the cell mixture and RNA mixture datasets.

## Scripts for reproducing a broader methods comparison

The [script] folder contains scripts that can reproduce the analysis and figures from our preprint: [scRNA-seq mixology: towards better benchmarking of single cell RNA-seq protocols and analysis methods](https://www.biorxiv.org/content/early/2018/10/03/433102). 

***Note:*** The `ggtern` package, which has been used to generate the ternary plots, has known issues with recent versions of `ggplot` and the relevant code may be broken if you have updated the `ggplot` package.
