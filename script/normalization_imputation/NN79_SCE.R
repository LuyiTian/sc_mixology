library(scPipe)
library(scater)
sce79.1 <- create_sce_by_dir("/Volumes/singlecellRNAseq-1/NN79_16.11.17/SC1")
sce79.2 <- create_sce_by_dir("/Volumes/singlecellRNAseq-1/NN79_16.11.17/SC2")
sce79.3 <- create_sce_by_dir("/Volumes/singlecellRNAseq-1/NN79_16.11.17/SC3")
sce79.4 <- create_sce_by_dir("/Volumes/singlecellRNAseq-1/NN79_16.11.17/SC4")



unigenes <- union(rownames(sce79.1), rownames(sce79.2))
unigenes <- union(unigenes, rownames(sce79.3))
unigenes <- union(unigenes, rownames(sce79.4))
unicells <- union(colnames(sce79.1), colnames(sce79.2))
unicells <- union(unicells, colnames(sce79.3))
unicells <- union(unicells, colnames(sce79.4))

new_cnt <- matrix(0, nrow = length(unigenes), ncol = length(unicells), dimnames = list(unigenes, unicells))
new_cnt[rownames(sce79.1), colnames(sce79.1)] = counts(sce79.1)
new_cnt[rownames(sce79.2), colnames(sce79.2)] = new_cnt[rownames(sce79.2), colnames(sce79.2)]+ counts(sce79.2)
new_cnt[rownames(sce79.3), colnames(sce79.3)] = new_cnt[rownames(sce79.3), colnames(sce79.3)]+ counts(sce79.3)
new_cnt[rownames(sce79.4), colnames(sce79.4)] = new_cnt[rownames(sce79.4), colnames(sce79.4)]+ counts(sce79.4)

new_qc <- as.matrix(colData(sce79.4))
new_qc <- new_qc + as.matrix(colData(sce79.3)[rownames(new_qc),])
new_qc <- new_qc + as.matrix(colData(sce79.2)[rownames(new_qc),])
new_qc <- new_qc + as.matrix(colData(sce79.1)[rownames(new_qc),])
new_qc <- data.frame(new_qc)

sce79 <- SingleCellExperiment(assays = list(counts = new_cnt))
QC_metrics(sce79) <- new_qc

demultiplex_info(sce79) <- demultiplex_info(sce79.4)

isSpike(sce79, "ERCC") = grepl("^ERCC-", rownames(sce79))

#sce79 <- fix_colname(sce79)

sce79.1 <- calculate_QC_metrics(sce79.1)
sce79.2 <- calculate_QC_metrics(sce79.2)
sce79.3 <- calculate_QC_metrics(sce79.3)
sce79.4 <- calculate_QC_metrics(sce79.4)
sce79 <- calculate_QC_metrics(sce79)
NN79_mRNA_annotation <- read.csv("FACS_design_benchmark.csv", stringsAsFactors = FALSE)
NN79_mRNA_annotation = NN79_mRNA_annotation[match(colnames(sce79.1),NN79_mRNA_annotation$wall_position),]

colData(sce79.1) = cbind(colData(sce79.1), DataFrame(NN79_mRNA_annotation))
colData(sce79.2) = cbind(colData(sce79.2), DataFrame(NN79_mRNA_annotation))
colData(sce79.3) = cbind(colData(sce79.3), DataFrame(NN79_mRNA_annotation))
colData(sce79.4) = cbind(colData(sce79.4), DataFrame(NN79_mRNA_annotation))
colData(sce79) = cbind(colData(sce79), DataFrame(NN79_mRNA_annotation))

sce79.1_qc <- detect_outlier(sce79.1, type="low", comp = 2)
table(QC_metrics(sce79.1_qc)$outliers)
sce79.2_qc <- detect_outlier(sce79.2, type="low", comp = 2)
table(QC_metrics(sce79.2_qc)$outliers)
sce79.3_qc <- detect_outlier(sce79.3, type="low", comp = 2)
table(QC_metrics(sce79.3_qc)$outliers)
sce79.4_qc <- detect_outlier(sce79.4, type="low", comp = 2)
table(QC_metrics(sce79.4_qc)$outliers)
sce79_qc <- detect_outlier(sce79, type="low", comp = 2)
table(QC_metrics(sce79_qc)$outliers)
outlier <- data.frame(table(QC_metrics(sce79.1_qc)$outliers), table(QC_metrics(sce79.2_qc)$outliers), table(QC_metrics(sce79.3_qc)$outliers), table(QC_metrics(sce79.4_qc)$outliers), table(QC_metrics(sce79_qc)$outliers))
write.csv(outlier, file = "outlier_NN79.csv")

sce79.1_qc = remove_outliers(sce79.1_qc)
sce79.2_qc = remove_outliers(sce79.2_qc)
sce79.3_qc = remove_outliers(sce79.3_qc)
sce79.4_qc = remove_outliers(sce79.4_qc)
sce79_qc = remove_outliers(sce79_qc)

keep1 = (apply(counts(sce79.1_qc), 1, function(x) mean(x[x>0])) > 1)  # average count larger than 1
keep2 = (rowSums(counts(sce79.1_qc)>0) > 10)  
table(keep1&keep2)
sce79.1_qc = sce79.1_qc[(keep1 & keep2), ]

keep1 = (apply(counts(sce79.2_qc), 1, function(x) mean(x[x>0])) > 1)  # average count larger than 1
keep2 = (rowSums(counts(sce79.2_qc)>0) > 10)  
table(keep1&keep2)
sce79.2_qc = sce79.2_qc[(keep1 & keep2), ]

keep1 = (apply(counts(sce79.3_qc), 1, function(x) mean(x[x>0])) > 1)  # average count larger than 1
keep2 = (rowSums(counts(sce79.3_qc)>0) > 10)  
table(keep1&keep2)
sce79.3_qc = sce79.3_qc[(keep1 & keep2), ]

keep1 = (apply(counts(sce79.4_qc), 1, function(x) mean(x[x>0])) > 1)  # average count larger than 1
keep2 = (rowSums(counts(sce79.4_qc)>0) > 10)  
table(keep1&keep2)
sce79.4_qc = sce79.4_qc[(keep1 & keep2), ]

keep1 = (apply(counts(sce79_qc), 1, function(x) mean(x[x>0])) > 1.5)  # average count larger than 1
keep2 = (rowSums(counts(sce79_qc)>0) > 10)  
table(keep1&keep2)
sce79_qc = sce79_qc[(keep1 & keep2), ]

is.spike.sc1 <- grepl("^ERCC", rownames(sce79.1_qc))
is.spike.sc2 <- grepl("^ERCC", rownames(sce79.2_qc))
is.spike.sc3 <- grepl("^ERCC", rownames(sce79.3_qc))
is.spike.sc4 <- grepl("^ERCC", rownames(sce79.4_qc))

sce79.1_qc <- calculateQCMetrics(sce79.1_qc, feature_controls = list(ERCC = is.spike.sc1))
sce79.2_qc <- calculateQCMetrics(sce79.2_qc, feature_controls = list(ERCC = is.spike.sc2))
sce79.3_qc <- calculateQCMetrics(sce79.3_qc, feature_controls = list(ERCC = is.spike.sc3))
sce79.4_qc <- calculateQCMetrics(sce79.4_qc, feature_controls = list(ERCC = is.spike.sc4))

sce79.1_qc <- sce79.1_qc[,!is.na(sce79.1_qc$H1975)]
sce79.2_qc <- sce79.2_qc[,!is.na(sce79.2_qc$H1975)]
sce79.3_qc <- sce79.3_qc[,!is.na(sce79.3_qc$H1975)]
sce79.4_qc <- sce79.4_qc[,!is.na(sce79.4_qc$H1975)]
sce79_qc <- sce79_qc[,!is.na(sce79_qc$H1975)]

#remove unuseful variables

rm(new_cnt, new_qc, NN79_mRNA_annotation, outlier, is.spike.sc1, is.spike.sc2, is.spike.sc3, is.spike.sc4, keep1, keep2, unicells, unigenes)

save.image("SCE_NN79.RData")


