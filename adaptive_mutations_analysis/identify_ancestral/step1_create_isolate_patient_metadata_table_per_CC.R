
# This R script is used to create two-column files with isolate and patient Ids for each CC

metadata_file = "/project/directory/colonising_pairs_patient-matched_gwas_dataset.all_pairs.metadata.csv";

metadata_table = read.delim(metadata_file, sep = "\t", header = T)
dim(metadata_table)
# [1] 2627   22

# Keep only QCd isolates
metadata_table = subset(metadata_table, QC1 == "yes" & QC2 == "yes")
dim(metadata_table)
# [1] 2270   22

# Get list of unique CC=
clonal_complexes = unique(c(as.vector(metadata_table$CC1), as.vector(metadata_table$CC2)))

# For each CC, extract list of unique patients and their isolates

for(c in 1:length(clonal_complexes))
{
  cc_table = mat.or.vec(1,2)
  tmp = which(metadata_table$CC1 == clonal_complexes[c])
  length(tmp)
  patients_cc = unique(as.vector(metadata_table$patientID[tmp]))
  length(patients_cc)
  for(p in 1:length(patients_cc))
  {
    tmp = which(metadata_table$CC1 == clonal_complexes[c] & metadata_table$patientID == patients_cc[p])
    isolates_patient = unique(c(as.vector(metadata_table$isolateID1[tmp]), as.vector(metadata_table$isolateID2[tmp])))
    for(i in 1:length(isolates_patient))
    {
      cc_table = rbind(cc_table, c(patients_cc[p], isolates_patient[i]))
    }
  }
  output_file = paste("colonising_pairs_patient-matched_gwas_dataset.CC",clonal_complexes[c],".host_isolate_ids.csv",sep="")
  print(output_file)
  cc_table = cc_table[-1,]
  colnames(cc_table) = c("host_id","isolate_id")
  write.table(cc_table,file=output_file,col.names=T,row.names=F,sep="\t",quote=F)
}
