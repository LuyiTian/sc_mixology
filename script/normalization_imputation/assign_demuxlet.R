load("~/mixture/SCE_sc_newQC.RData")
# cel-seq
demuxlet <- read.delim("NN84_RPI4_demuxlet.best")
index_anno <- read.csv("/wehisan/general/user_managed/grpu_naik.s_1/NN45_18.05.17/barcode_annotation.csv")
index_anno_alter <- index_anno
# use N instead of C or G for the end of barcode
index_anno_alter$index <- sub("C$", "N", index_anno$index)
index_anno_alter$index <- sub("G$", "N", index_anno_alter$index)
demuxlet_match <- merge(index_anno_alter, demuxlet,  by.y = "BARCODE", by.x = "index", all.x = TRUE)
# use NN instead of CN or GN, combine the results, delete the rows with total read < 500
index_anno_alter$index <- sub("GN$", "NN", index_anno_alter$index)
index_anno_alter$index <- sub("CN$", "NN", index_anno_alter$index)
demuxlet_match2 <- merge(index_anno_alter, demuxlet,  by.y = "BARCODE", by.x = "index", all.x = TRUE)
# get the lines end with NN
demuxlet_match2 <- demuxlet_match2[grepl("NN$", demuxlet_match2$index),]
demuxlet_match <- rbind(demuxlet_match, demuxlet_match2)
demuxlet_match <- demuxlet_match[demuxlet_match$RD.TOTL >= 500,]
sce4_qc$cell_line <- demuxlet_match[match(colnames(sce4_qc), demuxlet_match$samplename ), "SNG.1ST"]

# drop-seq
demuxlet <- read.delim("NN86_demuxlet.best")
index_anno <- read.csv("/wehisan/general/user_managed/grpu_naik.s_1/NN86_02.02.18/N701/index_anno.csv", header = FALSE)
names(index_anno) <- c("samplename", "index", "nreads")
index_anno <- index_anno[index_anno$samplename %in% colnames(sce701_qc),]

demuxlet_match <- merge(index_anno, demuxlet,  by.y = "BARCODE", by.x = "index", all.x = FALSE)
sce701_qc$cell_line <- demuxlet_match[match(colnames(sce701_qc), demuxlet_match$samplename ), "SNG.1ST"]

#test
library(scran)
library(scater)
sce701_qc <- computeSumFactors(sce701_qc)
sce701_qc <- normalize(sce701_qc)
pdf("MDS/Drop-seq_demuxlet.pdf")
plotMDS(sce701_qc, colour_by = "cell_line")
dev.off()
