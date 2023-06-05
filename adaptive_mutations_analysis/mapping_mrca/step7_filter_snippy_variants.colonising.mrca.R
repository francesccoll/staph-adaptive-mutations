
# This R script is used to filter out Snippy variants found in all isolates of the same host, that is, artifacts of mapping to the host MRCA sequence


isolate_metadata_file = "/project/directory/data/colonising_pairs_patient-matched_gwas_dataset.all_isolates.metadata.csv";

ref = "nctc8325"
ref = "je2"
ref = "emrsa15"

snippy_variants = paste("/project/directory/snippy_mrca/colonising_pairs_patient-matched_gwas_dataset.mrca.snippy_snpEff_annotated_variants.",ref,".in_ref_VCF.txt",sep="");
snippy_variants_kept = paste("/project/directory/snippy_mrca/colonising_pairs_patient-matched_gwas_dataset.mrca.snippy_snpEff_annotated_variants.",ref,".in_ref_VCF.kept.txt",sep="");
snippy_variants_rm = paste("/project/directory/snippy_mrca/colonising_pairs_patient-matched_gwas_dataset.mrca.snippy_snpEff_annotated_variants.",ref,".in_ref_VCF.rm.txt",sep="");

metadata = read.csv(isolate_metadata_file, sep = "\t", header = T)
dim(metadata)
# [1] 3497    9

metadata = subset(metadata, QC == "yes")
dim(metadata)
# [1] 3060    9

patients_metadata = unique(as.vector(metadata$patientID))
length(patients_metadata)
# [1] 791

variants = read.delim(snippy_variants, sep = "\t", header = T)
dim(variants)
# [1] 164984     29 > nctc8325
# [1] 123418     29 > je2
# [1] 145871     29 > emrsa15

#### Mutations with Ns in the reference or alternative allele will be deleted
tmp = which(grepl("N",variants$ref)==TRUE | grepl("N",variants$alt)==TRUE)
if(length(tmp)>0){ variants = variants[-tmp,]; }
dim(variants)
# [1] 164984     29 > nctc8325
# [1] 123418     29 > je2

# For variants with multiple snpEff annotations, remove intragenic annotations and keep only one annotation per variant
variants_to_remove_idx = vector(); # vector to store row indices to remove
samples = unique(as.vector(variants$sample))
for(s in 1:length(samples))
{
  print(s)
  # get variants for sample s
  tmp_s = which(variants$sample == samples[s])
  # get variants positions for sample s
  pos_s = unique(as.vector(variants$pos[tmp_s]))
  # for each variant position
  for(p in 1:length(pos_s))
  {
    tmp_p = which(variants$sample == samples[s] & variants$pos == pos_s[p])
    # if variant has multiple annotations
    if(length(tmp_p) > 1)
    {
      # remove intragenic_variant annotations
      tmp_p_rm = which(variants$sample == samples[s] & variants$pos == pos_s[p] & variants$annotation == "intragenic_variant")
      # make sure only one annotation is left after removing intragenic_variant annotations
      # NOTE: the same variant can lead to mutations in multiple locus_tag (if overlapping CDS)
      if((length(tmp_p) - length(tmp_p_rm))!=1)
      {
        tmp_p_keep = which(variants$sample == samples[s] & variants$pos == pos_s[p] & variants$annotation != "intragenic_variant")
        locus_tags = unique(as.vector(variants$locus_tag[tmp_p_keep]))
        if(length(locus_tags)==1)
        {
          print(paste("Error in position ",p," in sample ",s,sep=""));
        }
      }
      variants_to_remove_idx = c(variants_to_remove_idx, tmp_p_rm)
    }
  }
}

length(variants_to_remove_idx)
# [1] 99189
# [1] 55439

variants = variants[-variants_to_remove_idx,]
dim(variants)
# [1] 65795    29
# [1] 67979    29
# [1] 78996    29

get_patient_id = function(x){ y = paste(unlist(strsplit(unlist(strsplit(as.character(x),"[.]"))[-(1:3)] ,"-"))[1:2],collapse="-"); return(y); }
get_isolate_id = function(x){ y = paste(unlist(strsplit(unlist(strsplit(as.character(x),"[.]"))[1] ,"-"))[1],collapse="-"); return(y); }

variants$patient_id = sapply(variants$sample, get_patient_id)
variants$isolate_id = sapply(variants$sample, get_isolate_id)

patients_variants = unique(as.vector(variants$patient_id))
length(patients_variants)
# [1] 664
# [1] 705
# [1] 698

# check: make sure patients_variants are all in patients_metadata

length(which(is.na(match(patients_variants,patients_metadata))))
# [1] 0

# For each patient, extract isolate ids, and, for each variant, keep track of variants present in all patient's isolates (which will be filtered out)

variants_to_rm_idx = vector(); # vector to store variants to remove

for(p in 1:length(patients_variants))
{
  print(p)
  patient_id = as.character(patients_variants[p])
  patient_isolates = as.vector(metadata$isolateID[which(metadata$patientID==patient_id)])
  patient_variants_pos = unique(variants$pos[which(variants$patient_id == patient_id)])
  for(i in 1:length(patient_variants_pos))
  {
    variant_idx = which(variants$patient_id == patient_id & variants$pos == patient_variants_pos[i])
    variant_isolates = unique(as.vector(variants$isolate_id[variant_idx]))
    # if variant present in all patient's isolates, remove
    if(length(which(is.na(match(patient_isolates, variant_isolates))))==0)
    {
      variants_to_rm_idx = c(variants_to_rm_idx, variant_idx)
    }
  }
}
length(variants_to_rm_idx)
# [1] 8392
# [1] 8150 > je2
# [1] 9443 > emrsa15

variants = variants[-variants_to_rm_idx,]
dim(variants)
# [1] 57403    31
# [1] 59829    31 > je2
# [1] 69553    31 > emrsa15


write.table(variants, file=snippy_variants_kept, sep = '\t', col.names = T, row.names = F, quote = F)


