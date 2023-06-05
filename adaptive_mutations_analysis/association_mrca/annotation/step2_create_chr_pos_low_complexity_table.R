

dustmasker_output_file="/Users/francesccoll/fellowship/1.clinical_phenotype/8.patient-matched_col_inv/association/low_complexity/CC8_NCTC8325.dustmasker.txt"
chr_table_file="/Users/francesccoll/fellowship/1.clinical_phenotype/8.patient-matched_col_inv/association/cds/CC8_NCTC8325.cds_200bp_downstream.txt.gz"
output_file="/Users/francesccoll/fellowship/1.clinical_phenotype/8.patient-matched_col_inv/association/low_complexity/CC8_NCTC8325.dustmasker.chr_pos.txt"

dustmasker_output_file="/Users/francesccoll/fellowship/1.clinical_phenotype/8.patient-matched_col_inv/association/low_complexity/CC8_USA300_JE2.dustmasker.txt"
chr_table_file="/Users/francesccoll/fellowship/1.clinical_phenotype/8.patient-matched_col_inv/association/regulons_networks_ontologies/CC8_USA300_JE2.seif2019.TableS1.metabolic_submodules_subsystems.chr_pos.txt.gz"
output_file="/Users/francesccoll/fellowship/1.clinical_phenotype/8.patient-matched_col_inv/association/low_complexity/CC8_USA300_JE2.dustmasker.chr_pos.txt"


dustmasker_output_file="/Users/francesccoll/fellowship/1.clinical_phenotype/8.patient-matched_col_inv/association/low_complexity/CC22_EMRSA15.dustmasker.txt"
chr_table_file="/Users/francesccoll/fellowship/1.clinical_phenotype/8.patient-matched_col_inv/association/cds/CC22_EMRSA15.cds_200bp_downstream.txt.gz"
output_file="/Users/francesccoll/fellowship/1.clinical_phenotype/8.patient-matched_col_inv/association/low_complexity/CC22_EMRSA15.dustmasker.chr_pos.txt"


lmD = read.delim(dustmasker_output_file, sep = " ", header = T)
dim(lmD)
# [1] 1948    3, CC8_NCTC8325.dustmasker.txt
# [1] 2059    3, CC8_USA300_JE2.dustmasker.txt
# [1] 2002    3, CC22_EMRSA15.dustmasker.txt

window_lc = 5; # window around low complexity region to label as low complexity region as well

low_complexity_positions = vector()
for(r in 1:nrow(lmD))
{
  from = as.numeric(lmD[r,1]); to = as.numeric(lmD[r,3]); from = from - window_lc; to = to + window_lc;
  low_complexity_positions = c(low_complexity_positions, seq(from, to, 1))
}
low_complexity_positions = unique(sort(low_complexity_positions))
length(low_complexity_positions)
# [1] 69470
# [1] 74574
# [1] 71932



chr_table = read.delim(gzfile(chr_table_file) , sep = "\t", header = T)
chr_table = chr_table[,c(1:2)]
dim(chr_table)
# [1] 2821360       2, CC8_NCTC8325
# [1] 2874399       2, CC8_USA300_JE2
# [1] 2832298       2, CC22_EMRSA15

# Adding low complexity position to chromosome table
low_complexity = rep("no", nrow(chr_table))
low_complexity[low_complexity_positions] = "yes"
chr_table = cbind(chr_table, low_complexity)

write.table(chr_table, file=output_file, sep = '\t', col.names = T, row.names = F, quote = F)





