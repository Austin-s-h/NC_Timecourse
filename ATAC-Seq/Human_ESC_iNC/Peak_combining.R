library(readr)
library(BSgenome.Hsapiens.UCSC.hg38)
ESC_SampleSheet_Master <- read_csv("ATAC-Seq/Human_ESC_iNC/ESC_SampleSheet_Master.csv")
p1 <- import(ESC_SampleSheet_Master$Peaks[1])
p2 <- import(ESC_SampleSheet_Master$Peaks[2])
p3 <- import(ESC_SampleSheet_Master$Peaks[3])
p4 <- import(ESC_SampleSheet_Master$Peaks[4])
p5 <- import(ESC_SampleSheet_Master$Peaks[5])
p6 <- import(ESC_SampleSheet_Master$Peaks[6])
p7 <- import(ESC_SampleSheet_Master$Peaks[7])
p8 <- import(ESC_SampleSheet_Master$Peaks[8])
p9 <- import(ESC_SampleSheet_Master$Peaks[9])
p10 <- import(ESC_SampleSheet_Master$Peaks[10])
p11 <- import(ESC_SampleSheet_Master$Peaks[11])
p12 <- import(ESC_SampleSheet_Master$Peaks[12])

grangel <- GRangesList(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
hg38_seqinfo <- BSgenome.Hsapiens.UCSC.hg38@seqinfo
seqlevels(grangel) <- seqlevels(hg38_seqinfo)
seqinfo(grangel) <- hg38_seqinfo
combined <- c(unlist(grangel))
combined <- unique(combined)
combined <- resize(combined, 400, fix="center", use.names=TRUE,
             ignore.strand=T)
drop_ranges <-c("chrUn_KI270442v1","chrUn_KI270336v1","chrUn_KI270751v1",
"chr14_GL000225v1_random","chr22_KI270731v1_random", "chrUn_KI270519v1","chrUn_KI270315v1")
combined <- dropSeqlevels(combined, drop_ranges, pruning.mode = "coarse")
combined <- trim(combined)

export.bed(combined, con = "ATAC-Seq/Human_ESC_iNC/400bp_centered_all_peaks.bed")


