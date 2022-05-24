library(readr)
library(qdapTools)
output_global_TFs <- read_csv("~/workdir/diffTF/ATAC_Timecourse/timecourseDesign/manual_analysis/output.global.TFs.csv")
output_global_TFs_orig <- read_csv("~/workdir/diffTF/ATAC_Timecourse/timecourseDesign/manual_analysis/output.global.TFs.orig.csv")
all_timepoints_summary <- read_delim("~/workdir/diffTF/ATAC_Timecourse/timecourseDesign/FINAL_OUTPUT/extension100/all.timepoints.summary.tsv", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)
trans_table <- read_delim(file = "~/workdir/diffTF/src/TF_Gene_TranslationTables/JASPAR2020/translationTable_gg6.csv", delim = " ")
# https://www.grnpedia.org/trrust/
trrust_rawdata_human <- read_delim("~/workdir/diffTF/ATAC_Timecourse/timecourseDesign/manual_analysis/trrust_rawdata.human.tsv", 
                                   "\t", escape_double = FALSE, col_names = FALSE, 
                                   trim_ws = TRUE)
ddsTC_Condition_WE_vc_NC_All <- read_csv("~/local_git/NC_Timecourse/RNA-Seq/Exports/ddsTC_Condition_WE_vc_NC_All.csv")

all_timepoints_summary$Symbol <- lookup(all_timepoints_summary$TF, as.data.frame(trans_table[,c(3,1)]))
output_global_TFs$Symbol <- lookup(output_global_TFs$TF, as.data.frame(trans_table[,c(3,1)]))
output_global_TFs_orig$Symbol <- lookup(output_global_TFs_orig$TF, as.data.frame(trans_table[,c(3,1)]))

summary(as.factor(all_timepoints_summary$classification_q0.05_final))

activators <- all_timepoints_summary[all_timepoints_summary$classification_q0.05_final == "activator",]
repressors <- all_timepoints_summary[all_timepoints_summary$classification_q0.05_final == "repressor",]


represented <- all_timepoints_summary$Symbol[all_timepoints_summary$Symbol %in% trrust_rawdata_human$X1]
not_represented <- all_timepoints_summary$Symbol[!all_timepoints_summary$Symbol %in% trrust_rawdata_human$X1]

report <- data.frame()

for (i in seq(1,length(represented))) {
  db_sub <- subset(trrust_rawdata_human, trrust_rawdata_human$X1 == represented[i])
  act_count <- length(db_sub[db_sub$X3 == "Activation",4]$X4)
  rep_count <- length(db_sub[db_sub$X3 == "Repression",4]$X4)
  unkn_count <- length(db_sub[db_sub$X3 == "Unknown",4]$X4)
  report <- rbind(report, data.frame(represented[i], act_count, rep_count, unkn_count))
}

timepoints_by_symb <- as.data.frame(all_timepoints_summary[,c("Symbol","classification_q0.05_final")])
timepoints_by_symb <- timepoints_by_symb[!is.na(timepoints_by_symb$Symbol),]
timepoints_by_symb <- timepoints_by_symb[!duplicated(timepoints_by_symb$Symbol),]
report$NC_timecourse_call <- qdapTools::lookup(report$represented.i., timepoints_by_symb)
report$conflict <- ifelse(report$NC_timecourse_call == "activator" & report$rep_count > report$act_count, 
                          "conflict", "expected")

