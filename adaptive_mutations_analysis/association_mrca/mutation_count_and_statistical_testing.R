
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

# Mutations table
filtered_mutations_table = paste(dir,"mapping_mrca/",dataset,".mrca.snippy_snpEff_annotated_variants.",ref,".in_ref_VCF.filtered.txt.gz", sep = "");

# Metadata table
metadata_table = paste(dir,"data/",dataset,".all_isolates.metadata.csv", sep = "");

# Annotation files
gff_table_file = paste(dir,"mapping_mrca/",ref_file,".gff.table.txt", sep = "");
mader2016_ann_table_file = paste(dir,"association_mrca/annotation/",ref_file,".mader2016.annotation.txt.gz", sep = "");
mader2016_tu_table_file = paste(dir,"association_mrca/annotation/",ref_file,".operon_table.mader2016.txt", sep = "");
mader2016_tu_cluster_table_file = paste(dir,"association_mrca/annotation/",ref_file,".operon_cluster_table.mader2016.txt", sep = "");
mader2016_regulons_table_file = paste(dir,"association_mrca/annotation/",ref_file,".regulon_table.mader2016.txt", sep = "");
submodule_ann_file = paste(dir,"association_mrca/annotation/",ref_file,".seif2019.TableS1.metabolic_submodules.csv", sep = "");

# Change unit of association
feature_type = "CDS_only_functional"; # only mutations with HIGH or MODERATE functional annotation at CDS are included
feature_type = "TU_promoter"; # putative TU promoter region 200 bp downstream
feature_type = "TU_CDS_only_functional"; # only mutations with HIGH or MODERATE functional annotation at CDS are included
feature_type = "TU_cluster_CDS_only_functional"; # only mutations with HIGH or MODERATE functional annotation at CDS are included
feature_type = "TF_regulon_CDS_only_functional"; # only mutations with HIGH or MODERATE functional annotation at CDS are included
feature_type = "submodule_CDS_only_functional"; # metabolic submodules as defined by seif2019, defined only on JE2 reference genome

# Other parameteres
min_UA_length = 300; # minimum length of unit of association (in base pairs)
mutations_count_table = paste(dataset,".",feature_type,".",ref,".in_ref_VCF.mutations_count.txt",sep=""); 
mutations_count_table_assoc = paste(dataset,".",feature_type,".",ref,".in_ref_VCF.mutations_count.assoc.txt",sep="");



###############################################################################################################
####                                              FUNCTIONS                                                ####
###############################################################################################################

# Taken from https://github.com/johnlees/paired-samples/blob/master/poisson_test.R

poisson_test <- function(mutations, length, total_mutations, genome_length)
{
  print(length)
  length <- as.numeric(length)
  mutations <- as.numeric(mutations)
  p.value<-poisson.test(mutations,r=total_mutations*(length/genome_length), alternative = "greater")['p.value']
  return(as.numeric(p.value))
}


###############################################################################################################
####                                              MAIN CODE                                                ####
###############################################################################################################

mutations = read.delim(gzfile(filtered_mutations_table), sep = "\t")
dim(mutations)
# [1] 5552   32 > nctc8325, in_ref_VCF
# [1] 6257   32 > je2, in_ref_VCF
# [1] 6376   32 > emrsa15, in_ref_VCF

metadata = read.csv(metadata_table, sep = "\t", header = T)
dim(metadata)
# [1] 3497    9

#### If feature_type contains "CDS_only", then keep only mutations on CDS (intergenic are removed)
if(grepl("CDS_only",feature_type))
{
  mutations = subset(mutations, feature_type=="transcript")
}
dim(mutations)
# [1] 4550   32 > nctc8325, in_ref_VCF, CDS_only_functional
# [1] 5183   32 > je2, in_ref_VCF, CDS_only_functional
# [1] 5137   32 >  emrsa15, in_ref_VCF, CDS_only_functional

#### If feature_type contains "CDS_only", then keep only mutations on CDS (intergenic are removed)

if(grepl("_functional",feature_type))
{
  mutations = subset(mutations, annotation_impact=="HIGH" | annotation_impact=="MODERATE")
}
dim(mutations)
# [1] 3244   32 > nctc8325, in_ref_VCF, CDS_only_functional
# [1] 3549   32 > je2, in_ref_VCF, CDS_only_functional
# [1] 3446   32 > emrsa15, in_ref_VCF, CDS_only_functional

if(grepl("_promoter",feature_type))
{
  mutations = subset(mutations, annotation_impact=="MODIFIER")
}
dim(mutations)
# [1] 1006   32 > nctc8325, in_ref_VCF, TU_promoter

##### Saving mutations at locus of interest

locus_id = "SAOUHSC_02564";
mutationsL = mutations[which(mutations$locus_tag == locus_id),]
dim(mutationsL)
locus_output_file = paste(dataset,".",feature_type,".",ref,".",locus_id,".mutations.txt",sep="");
write.table(mutationsL, file=locus_output_file, sep = '\t', col.names = T, row.names = F, quote = F)

if(feature_type=="TU_promoter")
{
  locus_id = "U1094_promoter_200bp";
  mutationsL = mutations[which(grepl(locus_id, mutations$TU_promoter_Id_mader2016) == T),]
  dim(mutationsL)
  locus_output_file = paste(dataset,".",feature_type,".",ref,".",locus_id,".mutations.txt",sep="");
  write.table(mutationsL, file=locus_output_file, sep = '\t', col.names = T, row.names = F, quote = F)
}




##### Creating a table with mutations counts per loci

## Loading CDS annotation file
if(grepl("^submodule_CDS_only",feature_type))
{
  cds_table_file = paste(dir,"association_mrca/annotation/CC8_USA300_JE2.seif2019.TableS1.metabolic_submodules_subsystems.chr_pos.txt.gz", sep = ""); 
  cds_table = read.delim(gzfile(cds_table_file) , sep = "\t", header = T)
  dim(cds_table)
  # [1] 2874399       4
}

gene_table = read.delim(gff_table_file, sep = "\t", header = T)
dim(gene_table)
# [1] 2892    8, nctc8325
# [1] 2938    8, je2

mader2016_ann_table = read.delim(gzfile(mader2016_ann_table_file) , sep = "\t", header = T)
dim(mader2016_ann_table)
# [1] 2821360      19

## Getting list of unique loci Ids depending on unit of association selected
if(grepl("^CDS_only",feature_type))
{
  loci = unique(as.vector(gene_table$locus_tag))
  length(loci)
  # [1] 2892, nctc8325
  # [1] 2938, je2
}
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type))
{
  loci = unique(unlist(strsplit(as.vector(mader2016_ann_table$TU_Id_mader2016),";"))); # number of unique TU
  tmp = which(loci=="-"); if(length(tmp)>0){ loci = loci[-tmp]; }
  length(loci)
  # [1] 1418, nctc8325
}
if(grepl("^TU_promoter",feature_type))
{
  loci = unique(unlist(strsplit(as.vector(mader2016_ann_table$TU_promoter_Id_mader2016),";"))); # number of unique TU
  tmp = which(loci=="-"); if(length(tmp)>0){ loci = loci[-tmp]; }
  length(loci)
  # [1] 1418, nctc8325
}
if(grepl("^TU_cluster",feature_type))
{
  loci = unique(unlist(strsplit(as.vector(mader2016_ann_table$ClcortreeUA_mader2016),";"))); # number of unique TU clusters
  tmp = which(loci=="-"); if(length(tmp)>0){ loci = loci[-tmp]; }
  length(loci)
  # [1] 322, nctc8325
}
if(grepl("^TF_regulon",feature_type))
{
  loci = unique(unlist(strsplit(as.vector(mader2016_ann_table$transcription_factor_mader2016),";"))); # number of unique TF regulons
  tmp = which(loci=="-" | loci=="" | loci=="NA" | is.na(loci)); if(length(tmp)>0){ loci = loci[-tmp]; }
  length(loci)
  # [1] 47, nctc8325
}
if(grepl("^submodule_CDS_only",feature_type))
{
  loci = unique(unlist(strsplit(as.vector(cds_table$submodule),";"))); # number of unique CDS
  tmp = which(loci=="-"); if(length(tmp)>0){ loci = loci[-tmp]; }
  length(loci)
  # [1] 61
}



## Adding locus_id to mutations table
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type))
{
  xxx = match(mutations$pos, mader2016_ann_table$pos);
  mutations$TU_Id_mader2016 = mader2016_ann_table$TU_Id_mader2016[xxx]
  mutations$TU_Id_mader2016 = paste(mutations$TU_Id_mader2016,";",sep="")
}
if(grepl("^TU_promoter",feature_type))
{
  xxx = match(mutations$pos, mader2016_ann_table$pos);
  mutations$TU_promoter_Id_mader2016 = mader2016_ann_table$TU_promoter_Id_mader2016[xxx]
  mutations$TU_promoter_Id_mader2016 = paste(mutations$TU_promoter_Id_mader2016,";",sep="")
}
if(grepl("^TU_cluster",feature_type))
{
  xxx = match(mutations$pos, mader2016_ann_table$pos);
  mutations$ClcortreeUA_mader2016 = mader2016_ann_table$ClcortreeUA_mader2016[xxx]
  mutations$ClcortreeUA_mader2016 = paste(mutations$ClcortreeUA_mader2016,";",sep="")
}
if(grepl("^TF_regulon",feature_type))
{
  xxx = match(mutations$pos, mader2016_ann_table$pos);
  mutations$transcription_factor_mader2016 = mader2016_ann_table$transcription_factor_mader2016[xxx]
}
if(grepl("^submodule_CDS_only",feature_type))
{
  xxx = match(mutations$pos, cds_table$pos);
  mutations$locus_id = cds_table$submodule[xxx]
}
if(ref_file=="CC8_USA300_JE2")
{
  mutations$locus_tag = mutations$gene_name
}


## Saving mutations at locus of interest

if(feature_type=="TU_CDS_only_functional")
{
  locus_id = "U942";
  mutationsL = mutations[which(grepl(locus_id, mutations$TU_Id_mader2016) == T),]
  dim(mutationsL)
  locus_output_file = paste(dataset,".",feature_type,".",ref,".",locus_id,".mutations.txt",sep="");
  write.table(mutationsL, file=locus_output_file, sep = '\t', col.names = T, row.names = F, quote = F)
}


## Counting mutations per loci
loci_count_table = mat.or.vec(length(loci), 8); # table to keep locus name, mutations count, locus length, samples count
colnames(loci_count_table) = c('locus_id','feature_type','mutations_count','patient_count','patient_id','samples_count','samples_id','isolate_id');

for(c in 1:length(loci))
{
  print(c)
  locus_id = as.character(loci[c]); locus_id_grep = paste(locus_id,";",sep="");
  if(feature_type=="CDS"){ mutations_idx = which(grepl(locus_id_grep,mutations$locus_id)); }
  if(grepl("^CDS_only",feature_type)){ mutations_idx = which(mutations$locus_tag == locus_id); }
  if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type)){ mutations_idx = which(grepl(locus_id_grep,mutations$TU_Id_mader2016)); }
  if(grepl("^TU_promoter",feature_type)){ mutations_idx = which(grepl(locus_id_grep,mutations$TU_promoter_Id_mader2016)); }
  if(grepl("^TU_cluster",feature_type)){ mutations_idx = which(grepl(locus_id_grep,mutations$ClcortreeUA_mader2016)); }
  if(grepl("^TF_regulon",feature_type)){ mutations_idx = which(grepl(locus_id_grep,mutations$transcription_factor_mader2016)); }
  if(grepl("slice",feature_type)){ mutations_idx = which(grepl(locus_id_grep,mutations$locus_id)); }
  if(grepl("^submodule_CDS_only",feature_type)){ mutations_idx = which(grepl(locus_id_grep,mutations$locus_id)); }
  mutations_count = length(mutations_idx)
  samples_count = length(unique(mutations[mutations_idx,"sample"])); samples_id = paste(unique(mutations[mutations_idx,"sample"]),collapse = ";");
  patient_count = length(unique(mutations[mutations_idx,"patient_id"])); patient_id = paste(unique(mutations[mutations_idx,"patient_id"]),collapse = ";");
  isolate_id = paste(unique(mutations[mutations_idx,"isolate_id"]),collapse = ";");
  loci_count_table[c,] = c(locus_id, feature_type, mutations_count, patient_count, patient_id, samples_count, samples_id, isolate_id);
}

dim(loci_count_table)
# [1] 2892    8, nctc8325
# [1] 2938    8, je2
# [1] 2594    8, emrsa15
# [1] 1418    8, TU, nctc8325
# [1] 322   8, TU_cluster, nctc8325
# [1] 47  8, TF_regulon, nctc8325
# [1] 61  8, submodule_CDS_only, je2

write.table(loci_count_table, file = mutations_count_table, sep = '\t', col.names = T, row.names = F, quote = F);



###################################################################################################################
########                                          ADD LOCI ANNOTATION                                          ####
###################################################################################################################

if(grepl("^CDS_only",feature_type)){ loci_ann_table = read.delim(gff_table_file, sep = "\t", header = T); }
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type)){ loci_ann_table = read.delim(mader2016_tu_table_file, sep = "\t", header = T); }
if(grepl("^TU_cluster",feature_type)){ loci_ann_table = read.delim(mader2016_tu_cluster_table_file, sep = "\t", header = T); }
if(grepl("^TF_regulon",feature_type)){ loci_ann_table = read.delim(mader2016_regulons_table_file, sep = "\t", header = T); }
if(grepl("^submodule_CDS_only",feature_type)){ loci_ann_table = read.delim(submodule_ann_file, sep = "\t", header = T); }

dim(loci_ann_table)
# [1] 2892    8, nctc8325
# [1] 2938    8, je2
# [1] 1418   20, TU, nctc8325
# [1] 322  13, TU_cluster, nctc8325
# [1] 47  7, TF_regulon, nctc8325
# [1] 63 10, submodule_CDS_only, je2

loci_count_table = read.delim(mutations_count_table, sep = '\t', header = T)
dim(loci_count_table)
# [1] 2892    7, nctc8325
# [1] 2938    7, je2
# [1] 1418    7, TU, nctc8325
# [1] 322   7, TU_cluster, nctc8325
# [1] 61  8, submodule_CDS_only, je2

if(grepl("^CDS_only",feature_type)){ loci_count_table = merge(loci_count_table, loci_ann_table, by.x = "locus_id", by.y = "locus_tag"); }
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type)){ loci_count_table = merge(loci_count_table, loci_ann_table, by.x = "locus_id", by.y = "TU_Id_mader2016"); }
if(grepl("^TU_promoter",feature_type)){ loci_count_table = merge(loci_count_table, loci_ann_table, by.x = "locus_id", by.y = "TU_promoter_Id_mader2016"); }
if(grepl("^TU_cluster",feature_type)){ loci_count_table = merge(loci_count_table, loci_ann_table, by.x = "locus_id", by.y = "ClcortreeUA_mader2016"); }
if(grepl("^TF_regulon",feature_type)){ loci_count_table = merge(loci_count_table, loci_ann_table, by.x = "locus_id", by.y = "transcription_factor_mader2016"); }
if(grepl("^submodule_CDS_only",feature_type)){ loci_count_table = merge(loci_count_table, loci_ann_table, by.x = "locus_id", by.y = "submodule"); }

dim(loci_count_table)
# [1] 2892   14, nctc8325
# [1] 2938   14, je2
# [1] 1418   27, TU, nctc8325
# [1] 322  20, TU_cluster, nctc8325
# [1] 47 14, TU_cluster, nctc8325
# [1] 61 16, submodule_CDS_only, je2

if(grepl("^CDS_only",feature_type)){ loci_count_table$gene_length = (loci_count_table$end - loci_count_table$start)+1; }

###################################################################################################################
########                                          ASSOCIATION ANALYSIS                                         ####
###################################################################################################################

total_number_mutations = nrow(mutations)
if(ref == "nctc8325"){ chromosome_length = 2821360; }
if(ref == "je2"){ chromosome_length = 2874399; }
if(ref == "emrsa15"){ chromosome_length = 2832299; }
chromosome_num_loci = nrow(loci_count_table)


#### Removing loci with length shorter than min_UA_length (expect for TU_promoter, which are all 200bp)
if(grepl("^CDS_only",feature_type))
{
  loci_count_table = loci_count_table[-which(loci_count_table[,"gene_length"]<min_UA_length),]
}
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type) & !grepl("CDS_only",feature_type))
{
  loci_count_table = loci_count_table[-which(loci_count_table[,"TU_length_mader2016"]<min_UA_length),]
}
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type) & grepl("CDS_only",feature_type))
{
  loci_count_table = loci_count_table[-which(loci_count_table[,"TU_cds_length_mader2016"]<min_UA_length),]
}
if(grepl("^TU_cluster",feature_type) & !grepl("CDS_only",feature_type))
{
  loci_count_table = loci_count_table[-which(loci_count_table[,"ClcortreeUA_length_mader2016"]<min_UA_length),]
}
if(grepl("^TU_cluster",feature_type) & grepl("CDS_only",feature_type))
{
  loci_count_table = loci_count_table[-which(loci_count_table[,"ClcortreeUA_CDS_length_mader2016"]<min_UA_length),]
}
if(grepl("^TF_regulon",feature_type))
{
  tmp = which(loci_count_table[,"transcription_factor_length_mader2016"]<min_UA_length)
  if(length(tmp)>0){ loci_count_table = loci_count_table[-tmp,]; }
}
if(grepl("^submodule_CDS",feature_type))
{
  tmp = which(loci_count_table[,"submodule_cds_length"]<min_UA_length)
  if(length(tmp)>0){ loci_count_table = loci_count_table[-tmp,]; }
}



if(grepl("^CDS_only",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["gene_length"]], total_number_mutations, chromosome_length)})
}
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type) & !grepl("CDS_only",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["TU_length_mader2016"]], total_number_mutations, chromosome_length)})
}
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type) & grepl("CDS_only",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["TU_cds_length_mader2016"]], total_number_mutations, chromosome_length)})
}
if(grepl("^TU_promoter",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["TU_promoter_length_mader2016"]], total_number_mutations, chromosome_length)})
}
if(grepl("^TU_cluster",feature_type) & !grepl("CDS_only",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["ClcortreeUA_length_mader2016"]], total_number_mutations, chromosome_length)})
}
if(grepl("^TU_cluster",feature_type) & grepl("CDS_only",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["ClcortreeUA_CDS_length_mader2016"]], total_number_mutations, chromosome_length)})
}
if(grepl("^TF_regulon",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["transcription_factor_length_mader2016"]], total_number_mutations, chromosome_length)})
}
if(grepl("submodule_CDS_only",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["submodule_cds_length"]], total_number_mutations, chromosome_length)})
}

loci_count_table = cbind(loci_count_table, p_values)
# p_values_adj <- p.adjust(p_values, method="bonferroni", n=chromosome_num_loci)
p_values_adj <- p.adjust(p_values, method="BH", n=chromosome_num_loci)
loci_count_table = cbind(loci_count_table, p_values_adj)
p_value_cutoff = rep(as.character(0.05/chromosome_num_loci), nrow(loci_count_table))
loci_count_table = cbind(loci_count_table, p_value_cutoff)


# Add collection and CC samples came from
loci_count_table$collections = '-'
loci_count_table$CCs = '-'
for(r in 1:nrow(loci_count_table))
{
  samples = unlist(strsplit(x = as.character(loci_count_table$isolate_id[r]) , split = ";"))
  tmp = match(samples, metadata$isolateID)
  if(length(tmp)>0)
  {
    loci_count_table$collections[r] = as.character(paste(unique(metadata$collection[tmp]),collapse = ";"))
    loci_count_table$CCs[r] = as.character(paste(unique(metadata$CC[tmp]),collapse = ";"))
  }
}

write.table(loci_count_table, file=mutations_count_table_assoc, sep = '\t', col.names = T, row.names = F, quote = F)




