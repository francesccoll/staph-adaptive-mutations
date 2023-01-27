
## This R script is used to fit growth curves per replicate, to extract growth parameters (r and auc_l) and compare growth between isolates/strains

## Required R libraries
library(growthcurver)

## Reading and merging input data

rawdata = read.delim("curves_all_14-10-22.csv", sep = ",", header = T)
dim(rawdata)
# [1]  48 598
# Growthcurver vignette: The column containing the time must be named "time"
colnames(rawdata)[1] = "time"

# Growthcurver vignette: The remaining columns must have a unique well name that will be eventually be identified as the sample name
# NOTE: that means that replicates need to be averaged before fitting growth curves

metadata_curves = read.delim("curves_all_14-10-22.curves_metadata.csv", sep = "\t", header = T)
dim(metadata_curves)
# [1] 597   4

metadata_strains = read.delim("curves_all_14-10-22.strain_metadata.csv", sep = "\t", header = T)
dim(metadata_strains)
# [1] 35  4

metadata = merge(metadata_curves, metadata_strains, by.x = "sample_id", by.y = "strain_id")
dim(metadata)
# [1] 597   7
# adding label of mutation (e.g. vraA_p.Lys419fs_mutant)
# Important note: the variable label corresponds to sample_id but incorporates mutation annotation
metadata$sample_id = gsub("Absent", "wildtype", metadata$sample_id)
metadata$type = gsub("Absent", "wildtype", metadata$type)
metadata$label = paste(metadata$sample_id, metadata$gene, metadata$mutation, metadata$type, sep = "_")
metadata$mutation2 = paste(metadata$gene, metadata$mutation, sep = "_")

write.table(metadata, file = "curves_all_14-10-22.strain_curve_metadata.csv", sep = "\t", row.names = F, col.names = T, quote = F)


## For each curve_id (i.e. well), run growthcurver to extract growth parameters from fitted curve

# No daptomycin
environment = 0; # daptomycin concentration
from_time = 0; # >= 0 hours

curve_ids = unique(metadata$curve_id[which(metadata$environment == environment)])
curve_growth_par0 = mat.or.vec(length(curve_ids), 4)
colnames(curve_growth_par0) = c("curve_ids", "strain_id", "r", "auc_l")

for(c in 1:length(curve_ids)){
  print(c)
  sample_id = metadata$sample_id[which(metadata$curve_id == curve_ids[c])]
  subset = which(rawdata$time >= from_time)
  gc_fit <- SummarizeGrowth(rawdata$time[subset], rawdata[subset,curve_ids[c]])
  curve_growth_par0[c,] = c(curve_ids[c], sample_id, gc_fit$vals$r, gc_fit$vals$auc_l)
}

# With daptomycin
environment = 19
from_time = 7; # >= 7 hours

curve_ids = unique(metadata$curve_id[which(metadata$environment == environment)])
curve_growth_par19 = mat.or.vec(length(curve_ids), 4)
colnames(curve_growth_par19) = c("curve_ids", "strain_id", "r", "auc_l")

for(c in 1:length(curve_ids)){
  print(c)
  sample_id = metadata$sample_id[which(metadata$curve_id == curve_ids[c])]
  subset = which(rawdata$time >= from_time)
  gc_fit <- SummarizeGrowth(rawdata$time[subset], rawdata[subset,curve_ids[c]])
  curve_growth_par19[c,] = c(curve_ids[c], sample_id, gc_fit$vals$r, gc_fit$vals$auc_l)
}

write.table(curve_growth_par0, file = "curves_all_14-10-22.no_dap.parameters_fitted_curves.csv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(curve_growth_par19, file = "curves_all_14-10-22.with_dap.parameters_fitted_curves.csv", sep = "\t", row.names = F, col.names = T, quote = F)



## Comparison isolate growth parameters per mutation

mutations = unique(as.vector(metadata$mutation2))
mutations = mutations[-which(grepl("_knockout",mutations))]
mutations = mutations[-which(grepl("_none",mutations))]

environment = 0;
environment = 19;
if(environment == 0){
  curve_growth_par = curve_growth_par0; # choosing environment: daptomycin or not
}
if(environment == 19){
  curve_growth_par = curve_growth_par19; # choosing environment: daptomycin or not
}

results_table = mat.or.vec(1,12)
# NOTE: the p-values of one-way ANOVA and wilcox tests to compare means are saved, both for growth rate comparison (r) and AUC comparison (AUC)
colnames(results_table) = c("mutation","sample1","sample2","r1_mean","r2_mean","r_wt_p_value","r_aov_p_value","auc1_mean","auc2_mean","auc_wt_p_value","auc_aov_p_value","num_rep")

for(i in 1:length(mutations))
{
  mutation = mutations[i]
  print(paste(i, mutation, sep = " "))
  gene = unlist(strsplit(mutation, "_"))[1]
  # Extracting samples ids to be compared for mutation i
  samples_m = unique(metadata$sample_id[which(metadata$mutation2 == mutation & metadata$type == "mutant")]); # isolates with mutation
  samples_w = unique(metadata$sample_id[which(metadata$mutation2 == mutation & metadata$type == "wildtype")]); # wildtype isolates without mutation
  gene = unlist(strsplit(mutation, "_"))[1]
  sample_lockout = metadata$sample_id[which(metadata$gene == gene & metadata$type == "knockout")][1]; # gene transposon knockout
  sample_control = "JE2"; # transposon library strain - control
  
  # First, comparing mutant against wildtype isolates
  for(m in 1:length(samples_m)) {
    # getting growth paramters for mutant sample
    sample1 = samples_m[m]
    print(paste(m, sample1, sep = " "))
    r_m = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample1),"r"])
    r_m_mean = mean(r_m)
    auc_m = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample1),"auc_l"])
    auc_m_mean = mean(auc_m)
    num_rep = length(r_m)
    # getting growth rates for control/wildtype sample
    for(w in 1:length(samples_w)) {
      sample2 = samples_w[w]
      print(paste(w, sample2, sep = " "))
      r_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2),"r"])
      r_w_mean = mean(r_w)
      auc_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2),"auc_l"])
      auc_w_mean = mean(auc_w)
      # comparing growth parameters using wilcox.test
      test_r = wilcox.test(r_m, r_w, alternative = "greater", paired = FALSE, mu = 0)
      r_p_value_wt = test_r$p.value
      test_auc = wilcox.test(auc_m, auc_w, alternative = "greater", paired = FALSE, mu = 0)
      auc_p_value_wt = test_auc$p.value
      # comparing growth parameters using one-way ANOVA --> chosen as used in more papers, see growth_curves.analyses_and_plotting.docx, section 'Growth parameter comparisons performed in other papers'
      # preparing data format for ANOVA - r
      data = as.data.frame(cbind(c(r_m,r_w), c(rep("mutant",length(r_m)), rep("wildtype",length(r_m)))))
      colnames(data) = c("parameter","group")
      data$group <- ordered(data$group, levels = c("mutant", "wildtype"))
      data$parameter = as.numeric(data$parameter)
      aov <- aov(parameter ~ group, data = data)
      r_p_value_aov = summary(aov)[[1]][["Pr(>F)"]][1]
      # preparing data format for ANOVA - AUC
      data = as.data.frame(cbind(c(auc_m,auc_w), c(rep("mutant",length(auc_m)), rep("wildtype",length(auc_w)))))
      colnames(data) = c("parameter","group")
      data$group <- ordered(data$group, levels = c("mutant", "wildtype"))
      data$parameter = as.numeric(data$parameter)
      aov <- aov(parameter ~ group, data = data)
      auc_p_value_aov = summary(aov)[[1]][["Pr(>F)"]][1]
      # adding results as new row
      results_table = rbind(results_table, c(mutation, sample1, sample2, r_m_mean, r_w_mean, r_p_value_wt, r_p_value_aov, auc_m_mean, auc_w_mean, auc_p_value_wt, auc_p_value_aov, num_rep))
    }
    # Second, comparing mutant vs. transposon knockout
    sample2 = sample_lockout
    r_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2),"r"])
    r_w_mean = mean(r_w)
    auc_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2),"auc_l"])
    auc_w_mean = mean(auc_w)
    # comparing growth parameters using wilcox.test
    test_r = wilcox.test(r_m, r_w, alternative = "greater", paired = FALSE, mu = 0)
    r_p_value_wt = test_r$p.value
    test_auc = wilcox.test(auc_m, auc_w, alternative = "greater", paired = FALSE, mu = 0)
    auc_p_value_wt = test_auc$p.value
    # comparing growth parameters using one-way ANOVA --> chosen as used in more papers, see growth_curves.analyses_and_plotting.docx, section 'Growth parameter comparisons performed in other papers'
    # preparing data format for ANOVA - r
    data = as.data.frame(cbind(c(r_m,r_w), c(rep("mutant",length(r_m)), rep("wildtype",length(r_m)))))
    colnames(data) = c("parameter","group")
    data$group <- ordered(data$group, levels = c("mutant", "wildtype"))
    data$parameter = as.numeric(data$parameter)
    aov <- aov(parameter ~ group, data = data)
    r_p_value_aov = summary(aov)[[1]][["Pr(>F)"]][1]
    # preparing data format for ANOVA - AUC
    data = as.data.frame(cbind(c(auc_m,auc_w), c(rep("mutant",length(auc_m)), rep("wildtype",length(auc_w)))))
    colnames(data) = c("parameter","group")
    data$group <- ordered(data$group, levels = c("mutant", "wildtype"))
    data$parameter = as.numeric(data$parameter)
    aov <- aov(parameter ~ group, data = data)
    auc_p_value_aov = summary(aov)[[1]][["Pr(>F)"]][1]
    # adding results as new row
    results_table = rbind(results_table, c(mutation, sample1, sample2, r_m_mean, r_w_mean, r_p_value_wt, r_p_value_aov, auc_m_mean, auc_w_mean, auc_p_value_wt, auc_p_value_aov, num_rep))
    
    # Third, comparing mutant vs. control
    sample2 = sample_control
    r_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2),"r"])
    r_w_mean = mean(r_w)
    auc_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2),"auc_l"])
    auc_w_mean = mean(auc_w)
    # comparing growth parameters using wilcox.test
    test_r = wilcox.test(r_m, r_w, alternative = "greater", paired = FALSE, mu = 0)
    r_p_value_wt = test_r$p.value
    test_auc = wilcox.test(auc_m, auc_w, alternative = "greater", paired = FALSE, mu = 0)
    auc_p_value_wt = test_auc$p.value
    # comparing growth parameters using one-way ANOVA --> chosen as used in more papers, see growth_curves.analyses_and_plotting.docx, section 'Growth parameter comparisons performed in other papers'
    # preparing data format for ANOVA - r
    data = as.data.frame(cbind(c(r_m,r_w), c(rep("mutant",length(r_m)), rep("wildtype",length(r_m)))))
    colnames(data) = c("parameter","group")
    data$group <- ordered(data$group, levels = c("mutant", "wildtype"))
    data$parameter = as.numeric(data$parameter)
    aov <- aov(parameter ~ group, data = data)
    r_p_value_aov = summary(aov)[[1]][["Pr(>F)"]][1]
    # preparing data format for ANOVA - AUC
    data = as.data.frame(cbind(c(auc_m,auc_w), c(rep("mutant",length(auc_m)), rep("wildtype",length(auc_w)))))
    colnames(data) = c("parameter","group")
    data$group <- ordered(data$group, levels = c("mutant", "wildtype"))
    data$parameter = as.numeric(data$parameter)
    aov <- aov(parameter ~ group, data = data)
    auc_p_value_aov = summary(aov)[[1]][["Pr(>F)"]][1]
    # adding results as new row
    results_table = rbind(results_table, c(mutation, sample1, sample2, r_m_mean, r_w_mean, r_p_value_wt, r_p_value_aov, auc_m_mean, auc_w_mean, auc_p_value_wt, auc_p_value_aov, num_rep))
  }
}

if(environment == 0){
  output = "fitted_growth_curves_per_replicate.no_dap.growth_parameters.comparison.csv"
}
if(environment == 19){
  output = "fitted_growth_curves_per_replicate.with_dap.growth_parameters.comparison.csv"
}

write.table(results_table, file = output, sep = "\t", row.names = F, col.names = T, quote = F)




