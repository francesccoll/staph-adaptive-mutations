

###############################################################################################################
####                                         INPUT/OUTPUT FILES                                            ####
###############################################################################################################

dataset = "colonising_pairs_patient-matched_gwas_dataset";

# Change ref and ref_file variables to change reference genome
ref = "nctc8325"; ref_file = "CC8_NCTC8325";
ref = "je2"; ref_file = "CC8_USA300_JE2";
ref = "emrsa15"; ref_file = "CC22_EMRSA15";


setwd("/project/directory/association_mrca/");
dir = "/project/directory/"
combined_table = paste(dir,"mapping_mrca/",dataset,".mrca.snippy_snpEff_annotated_variants.",ref,".in_ref_VCF.kept.txt.gz", sep = "");
gff_table_file = paste(dir,"annotation/",ref_file,".gff.table.txt", sep = "");
low_complexity_file = paste(dir,"annotation/",ref_file,".dustmasker.chr_pos.txt.gz", sep = "");
repetitive_regions_file = paste(dir,"annotation/",ref_file,".repetitive_regions.tab", sep = "");

recombination_mutations_table = paste(dir,"mapping_mrca/",dataset,".mrca.snippy_snpEff_annotated_variants.",ref,".in_ref_VCF.recombination.txt", sep = "");
filtered_mutations_table = paste(dir,"mapping_mrca/",dataset,".mrca.snippy_snpEff_annotated_variants.",ref,".in_ref_VCF.filtered.txt", sep = "");

###############################################################################################################
####                                              MAIN CODE                                                ####
###############################################################################################################

mutations = read.delim(gzfile(combined_table), sep = "\t")
dim(mutations)
# [1] 57403    31 > nctc8325, in_ref_VCF
# [1] 59829    31 > je2, in_ref_VCF
# [1] 69553    31 > emrsa15, in_ref_VCF

###### Mutations at regions of low complexity need to be removed
lc_table = read.delim(gzfile(low_complexity_file) , sep = "\t", header = T)
dim(lc_table)
# [1] 2821360       3, nctc8325
# [1] 2874399       3, je2
# [1] 2832298       3, emrsa15

low_complexity = lc_table$low_complexity[mutations$pos]
mutations = cbind(mutations, low_complexity)

mutations = subset(mutations, low_complexity == "no")
dim(mutations)
# [1] 55268    32 > nctc8325, in_ref_VCF
# [1] 67564    32 > emrsa15, in_ref_VCF

###### Mutations at repetitive regions need to be removed
rep_table = read.table(repetitive_regions_file, sep = "\t", header = F)
dim(rep_table)
# [1] 409   2
rep_chr_pos = vector()
for(r in 1:nrow(rep_table))
{
  rep_chr_pos = c(rep_chr_pos, seq(rep_table[r,1], rep_table[r,2], 1))
}
length(unique(rep_chr_pos))
# [1] 212381 > nctc8325
# [1] 216470 > emrsa15
tmp = which(!is.na(match(mutations$pos, rep_chr_pos)))
if(length(tmp)>0){ mutations = mutations[-tmp,]; }
dim(mutations)
# [1] 51888    32 > nctc8325, in_ref_VCF
# [1] 63663    32 > emrsa15, in_ref_VCF


#### Hosts with a much higher number of mutations (outliers) than the rest are removed

outlier_hosts = c("price2016-H238")
mutations = mutations[-which(!is.na(match(mutations$patient_id, outlier_hosts))),]
dim(mutations)
# [1] 14047    32 > nctc8325, in_ref_VCF


#### Hosts with a much higher number of colonies/isolates sequenced than the rest will be removed 
outlier_hosts = c("tong2015-T126")
mutations = mutations[-which(!is.na(match(mutations$patient_id, outlier_hosts))),]
dim(mutations)
# [1] 13696    32 > nctc8325, in_ref_VCF


#### Number of unique samples (i.e. pairs) with mutations
laSamples = unique(as.vector(mutations$sample))
length(laSamples)
# [1] 1862 >  nctc8325, in_ref_VCF

#### Recombination events need to be removed

window = 1000; # a window of 100 failed to remove recombination in phage region
laRecomPos = vector(); # vector to store recombination positions to delete
for(s in 1:length(laSamples))
{
  print(s)
  sss = which(mutations$sample==as.character(laSamples[s]))
  mutationsSam = mutations[sss,];
  for(m in 1:length(sss))
  {
    recom_left_idx = which((mutationsSam$pos < mutationsSam$pos[m]) & (mutationsSam$pos >= (mutationsSam$pos[m]-window))); # check if neighbouring mutations downstream
    if(length(recom_left_idx)>0){ laRecomPos = c(laRecomPos, sss[recom_left_idx]); }
    recom_right_idx = which((mutationsSam$pos > mutationsSam$pos[m]) & (mutationsSam$pos <= (mutationsSam$pos[m]+window))); # check if neighbouring mutations upstream
    if(length(recom_right_idx)>0){ laRecomPos = c(laRecomPos, sss[recom_right_idx]); }
  }
}
laRecomPos = unique(laRecomPos);
length(laRecomPos)
# [1] 3922 > nctc8325, in_ref_VCF

mutationsRcb = mutations[laRecomPos,];
dim(mutationsRcb)
# [1] 3922   32 > nctc8325, in_ref_VCF

write.table(mutationsRcb, file=recombination_mutations_table, sep = '\t', col.names = T, row.names = F, quote = F)

if(length(laRecomPos)>0){ mutations = mutations[-laRecomPos,]; }
dim(mutations)
# [1] 9774   32 > nctc8325, in_ref_VCF


# mutations found in multiple isolates of the same patient need to be removed (de-duplicated)

get_patient_id = function(x){ y = paste(unlist(strsplit(unlist(strsplit(as.character(x),"[.]"))[-(1:3)] ,"-"))[1:2],collapse="-"); return(y); }
mutations$patient_id = sapply(mutations$sample, get_patient_id)
patients = unique(as.vector(mutations$patient_id))
length(patients)
# [1] 644 > nctc8325, in_ref_VCF


mutations_to_remove = vector();
for(p in 1:length(patients))
{
  tmp = which(mutations$patient_id == patients[p])
  tmp2 = which(duplicated(mutations$pos[tmp]))
  if(length(tmp2)>0){ mutations_to_remove = c(mutations_to_remove, tmp[tmp2]); }
}
mutations = mutations[-mutations_to_remove,];
dim(mutations)
# [1] 5856   32 > nctc8325, in_ref_VCF

# NOTE: a total of 5,856 mutations can be attributed to point mutations (as opposed to recombination)


### Removing mutations at phage regions

if(ref == "nctc8325")
{
  phage_regions_from = c(1923408, 1463618, 2031911)
  phage_regions_to = c(1967011, 1508580, 2074632)
}
if(ref == "je2")
{
  phage_regions_from = c(1563373, 2086236)
  phage_regions_to = c( 1592434, 2098004)
}
if(ref == "emrsa15")
{
  phage_regions_from = c(1533521, 2042935)
  phage_regions_to = c(1552285, 2054703)
}

phage_mutations_idx = vector()
for(p in 1:length(phage_regions_from))
{
  tmp = which(mutations$pos >= phage_regions_from[p] & mutations$pos <= phage_regions_to[p])
  phage_mutations_idx = c(phage_mutations_idx, tmp)
}

mutations = mutations[-phage_mutations_idx,]
dim(mutations)
# [1] 5552   32 > nctc8325, in_ref_VCF



write.table(mutations, file=filtered_mutations_table, sep = '\t', col.names = T, row.names = F, quote = F)






