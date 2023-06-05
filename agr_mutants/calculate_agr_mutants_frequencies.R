
# This R script is used to calculate agr-mutant frequencies in different subsets of patients

carriage_metadata_file = "/project/directory/data/colonising_pairs_patient-matched_gwas_dataset.all_pairs.metadata.xlsx"
carriage_mutations_file = "/project/directory/data/colonising_pairs_patient-matched_gwas_dataset.all_pairs.snippy_snpEff_annotated_variants.nctc8325.filtered.txt"

output = "colonising_patient.agr_summary.csv";

require(gdata)

carriage_metadata = read.xls(carriage_metadata_file, header = T, sheet = 1)
dim(carriage_metadata)
# [1] 2627   22


### keeping only isolates that passed QC

carriage_metadata = subset(carriage_metadata, QC1 == "yes" & QC2 == "yes")
dim(carriage_metadata)
# [1] 2275   22

carriage_hosts = unique(carriage_metadata$patientID)
length(carriage_hosts)
# [1] 793

# Extracting infected and non-infected patients

carriage_hosts_asymp = unique(carriage_metadata$patientID[which(carriage_metadata$Symptomatic=="no")])
length(carriage_hosts_asymp)
# [1] 497

carriage_hosts_symp = unique(carriage_metadata$patientID[which(carriage_metadata$Symptomatic=="yes")])
length(carriage_hosts_symp)
# [1] 206

carriage_hosts_unknown = unique(carriage_metadata$patientID[which(carriage_metadata$Symptomatic=="unknown")])
length(carriage_hosts_unknown)
# [1] 90

## Loading filtered mutations 

carriage_mutations = read.delim(carriage_mutations_file, sep = "\t", header = T)
dim(carriage_mutations)
# [1] 4691   31


### Preparing table with agr-mutant information per patient

agr_mutants_table = mat.or.vec(length(carriage_hosts),17)
colnames(agr_mutants_table) = c("patient","sites","collection","symptomatic","clonal_complex","isolates_num","annotation_impact_agrA","num_mutations_agrA","annotation_impact_agrC","num_mutations_agrC","agr_mutant","min_core_snps", "max_core_snps", "median_core_snps", "min_time_span", "max_time_span", "num_wVar")
dim(agr_mutants_table)
# [1] 793  17

for(p in 1:length(carriage_hosts))
{
  patient = as.character(carriage_hosts[p])
  tmp = which(carriage_metadata$patientID == patient)
  
  # Extracting agrA (SAOUHSC_02265) and agrC mutations (SAOUHSC_02264)
  annotation_impact_agrA = "NA"; num_mutations_agrA = "0"; 
  agrA_tmp = which(carriage_mutations$patient_id == patient & carriage_mutations$locus_tag == "SAOUHSC_02265")
  if(length(agrA_tmp)>0){ num_mutations_agrA = as.character(length(agrA_tmp)); annotation_impact_agrA = paste(unique(carriage_mutations$annotation_impact[agrA_tmp]), collapse = ";"); }
  annotation_impact_agrC = "NA"; num_mutations_agrC = "0"; 
  agrC_tmp = which(carriage_mutations$patient_id == patient & carriage_mutations$locus_tag == "SAOUHSC_02264")
  if(length(agrC_tmp)>0){ num_mutations_agrC = as.character(length(agrC_tmp)); annotation_impact_agrC = paste(unique(carriage_mutations$annotation_impact[agrC_tmp]), collapse = ";"); }
  # Determining agr-inactivating mutants
  agr_mutant = "no"
  if(grepl("HIGH", annotation_impact_agrA) | grepl("MODERATE", annotation_impact_agrA) | grepl("HIGH", annotation_impact_agrC) | grepl("MODERATE", annotation_impact_agrC)){ agr_mutant = "yes"; }
  
  collection = as.character(unique(carriage_metadata$collection[tmp]))
  symptomatic = as.character(unique(carriage_metadata$Symptomatic[tmp]))
  sites = paste(unique(c(as.vector(carriage_metadata$site1[tmp]),as.vector(carriage_metadata$site2[tmp]))), collapse = ";")
  clonal_complex = paste(unique(c(as.vector(carriage_metadata$CC1[tmp]),as.vector(carriage_metadata$CC2[tmp]))), collapse = ";")
  isolates = unique(c(as.vector(carriage_metadata$isolateID1[tmp]), as.vector(carriage_metadata$isolateID2[tmp])))
  isolates_text = paste(isolates, collapse = ";")
  isolates_num = as.character(length(isolates))
  core_snps = as.vector(carriage_metadata$SNPdisCoreGenome[tmp])
  min_core_snps = as.character(min(core_snps)); max_core_snps = as.character(max(core_snps)); median_core_snps = as.character(median(core_snps));
  # Time distance between isolates
  max_time_span = "not_calculated"; min_time_span = "not_calculated";
  if(collection=="paterson2015" | collection=="chow2017")
  {
    max_time_span = "unknown";
    min_time_span = "unknown";
  } else
  {
      time_distances = as.Date(carriage_metadata$collection_date2[tmp], format = "%Y-%m-%d") - as.Date(carriage_metadata$collection_date1[tmp], format = "%Y-%m-%d");
      max_time_span = as.character(max(as.integer(time_distances)));
      min_time_span = as.character(min(as.integer(time_distances)));
  }
  # Extracting number of whole-genome mutations
  tmp2 = which(carriage_mutations$patient_id == patient)
  num_wVar = as.character(length(tmp2))
  
  # Saving information
  newrow = c(patient, sites, collection, symptomatic, clonal_complex, isolates_num, annotation_impact_agrA, num_mutations_agrA, annotation_impact_agrC, num_mutations_agrC, agr_mutant, min_core_snps, max_core_snps, median_core_snps, min_time_span, max_time_span, num_wVar)
  agr_mutants_table[p,] = newrow
}


# Keeping only major CC
major_CCs = c("59", "398", "8", "15","5","1","45","30","22")
agr_mutants_table[-which(!is.na(match(agr_mutants_table[,"clonal_complex"],major_CCs))),"clonal_complex"] = "minor_ST"
dim(agr_mutants_table)


summary_table_file = "colonising_pairs_patient-matched_gwas_dataset.all_pairs.agr_summary.csv";
write.table(agr_mutants_table, file=summary_table_file, sep = '\t', col.names = T, row.names = F, quote = F)


##### Plotting the data




##### Logistic regression model to test the effect of several variables on the emergence of agr mutants

agr_mutants_table = read.delim(summary_table_file, sep = "\t", header = T)
dim(agr_mutants_table)
# [1] 793  17

# First remove patients of unknown symtomatic carrier status

agr_mutants_table = subset(agr_mutants_table, symptomatic != "unknown")
dim(agr_mutants_table)
# [1] 703  17

# Note: the min_core_snps was instead of using the max_core_snps as the latter is potentially 

mylogit <- glm(agr_mutant ~ symptomatic + isolates_num + min_core_snps + collection + clonal_complex, data = agr_mutants_table, family = "binomial")
summary(mylogit)

# Call:
#   glm(formula = agr_mutant ~ symptomatic + isolates_num + min_core_snps + 
#         collection + clonal_complex, family = "binomial", data = agr_mutants_table)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.5198  -0.3693  -0.2440  -0.1441   2.8635  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              -3.37560    1.31384  -2.569  0.01019 *  
# symptomaticyes            0.58953    0.42852   1.376  0.16890    
# isolates_num              0.18814    0.05004   3.759  0.00017 ***
# min_core_snps             0.05729    0.02365   2.423  0.01541 *  
# collectionharrison2016  -19.11327 2229.85003  -0.009  0.99316    
# collectionmpros          -0.40093    1.19580  -0.335  0.73741    
# collectionpaterson2015  -24.28897 3884.49742  -0.006  0.99501    
# collectionprice2013       0.05500    1.23573   0.045  0.96450    
# collectionprice2016      -1.42898    1.19599  -1.195  0.23216    
# collectionyoung2017      -2.14312    1.38895  -1.543  0.12284    
# clonal_complex15        -16.25753 1561.26566  -0.010  0.99169    
# clonal_complex22          0.52005    0.82522   0.630  0.52857    
# clonal_complex30         -0.27653    0.91307  -0.303  0.76200    
# clonal_complex398       -16.41096 2320.27315  -0.007  0.99436    
# clonal_complex45         -0.57394    1.26713  -0.453  0.65059    
# clonal_complex5          -0.54056    1.26826  -0.426  0.66995    
# clonal_complex59          0.29557    1.32581   0.223  0.82358    
# clonal_complex8           0.55796    1.06951   0.522  0.60188    
# clonal_complexminor_ST    0.22857    0.97952   0.233  0.81549    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 289.91  on 702  degrees of freedom
# Residual deviance: 250.05  on 684  degrees of freedom
# AIC: 288.05
# 
# Number of Fisher Scoring iterations: 18


# Calcularing odds rations and confidence intervals

exp(0.18814)
# [1] 1.207002
exp(0.18814 - 1.96 * 0.05004)
# [1] 1.094242
exp(0.18814 + 1.96 * 0.05004)
# [1] 1.331383

exp(confint(mylogit))
# Waiting for profiling to be done...
# 2.5 %       97.5 %
# (Intercept)             1.320388e-03 3.148109e-01
# symptomaticyes          7.591232e-01 4.135203e+00
# isolates_num            1.096670e+00 1.339259e+00
# min_core_snps           1.007364e+00 1.109006e+00

