
## This R script is used to fit growth curves per replicate, to extract growth parameters (r and auc_l) and compare growth between isolates/strains

## Required R libraries
library(growthcurver)

## Reading and merging input data

rawdata = read.delim("Export_All_Plates_Single_File.curves_growthdata.csv", sep = "\t", header = T)
dim(rawdata)
# [1]  193 2881

# Growthcurver vignette: The column containing the time must be named "time"
colnames(rawdata)[1]
# [1] "time"

# NOTE: curve_ids in rawdata (column names) need to be corrected
colnames(rawdata) = gsub("^X",  "", colnames(rawdata))

# Growthcurver vignette: The remaining columns must have a unique well name that will be eventually be identified as the sample name
# NOTE: that means that replicates need to be averaged before fitting growth curves

metadata_curves = read.delim("Export_All_Plates_Single_File.curves_metadata.csv", sep = "\t", header = T)
dim(metadata_curves)
# [1] 2880    6

metadata_strains = read.delim("Export_All_Plates_Single_File.strain_metadata.csv", sep = "\t", header = T)
dim(metadata_strains)
# [1] 10  6

metadata = merge(metadata_curves, metadata_strains, by.x = "sample_id", by.y = "strain_id")
dim(metadata)
# [1] 2880   11

# adding label of mutation (e.g. vraA_p.Lys419fs_mutant)
# Important note: the variable label corresponds to sample_id but incorporates mutation annotation

metadata$label = paste(metadata$isolate_id, metadata$gene, metadata$mutation, metadata$type, sep = "_")
metadata$mutation2 = paste(metadata$gene, metadata$mutation, sep = "_")
dim(metadata)
# [1] 2880   13

write.table(metadata, file = "Export_All_Plates_Single_File.strain_curve_metadata.csv", sep = "\t", row.names = F, col.names = T, quote = F)


## For each curve_id (i.e. well), run growthcurver to extract growth parameters from fitted curve

curve_ids = unique(metadata$curve_id)
length(curve_ids)
# [1] 2880
curve_growth_par = mat.or.vec(length(curve_ids), 5)
colnames(curve_growth_par) = c("curve_id", "strain_id", "nitrogen_source_pm3", "r", "auc_l")

for(c in 1:length(curve_ids)){
  print(c)
  sample_id = metadata$sample_id[which(metadata$curve_id == curve_ids[c])]
  nitrogen_source = metadata$nitrogen_source_pm3[which(metadata$curve_id == curve_ids[c])]
  #subset = which(rawdata$time >= from_time)
  #gc_fit <- SummarizeGrowth(rawdata$time[subset], rawdata[subset,curve_ids[c]])
  gc_fit <- SummarizeGrowth(rawdata$time, rawdata[,curve_ids[c]])
  curve_growth_par[c,] = c(curve_ids[c], sample_id, nitrogen_source, gc_fit$vals$r, gc_fit$vals$auc_l)
}

write.table(curve_growth_par, file = "Export_All_Plates_Single_File.parameters_fitted_curves.csv", sep = "\t", row.names = F, col.names = T, quote = F)

## Editing JE2 name

metadata$mutation2 = gsub("none_none", "JE2_control", metadata$mutation2)
metadata$label = gsub("JE2_none_none_knockout", "JE2_control", metadata$label)


## Comparison isolate growth parameters per mutation across nitrogen sources

results_table = mat.or.vec(1,13)
# NOTE: the p-values of one-way ANOVA and wilcox tests to compare means are saved, both for growth rate comparison (r) and AUC comparison (AUC)
colnames(results_table) = c("mutation","nitrogen_source","sample1","sample2","r1_mean","r2_mean","r_wt_p_value","r_aov_p_value","auc1_mean","auc2_mean","auc_wt_p_value","auc_aov_p_value","num_rep")

# 1. comparing fitted growth parameters of nasD and ureG transposons against controls (JE2 and spa) across all nitrogen sources

knockouts = c("nasD_knockout", "ureG_knockout")
knockouts_strain_ids = c("NE857", "NE1113")
controls = c("spa_knockout", "JE2_control")
controls_strain_ids = c("NE286", "JE2")

nitrogen_sources = unique(sort(metadata$nitrogen_source_pm3))

for(s in 1:length(nitrogen_sources))
{
  nitrogen_source = nitrogen_sources[s]
  print(paste(s, nitrogen_source, sep = " "))
  
for(i in 1:length(knockouts))
{
  knockout = knockouts[i]
  print(paste(i, knockout, sep = " "))
  
    # getting growth paramters for knockout sample
    sample1 = knockouts_strain_ids[i]
    print(paste(i, sample1, sep = " "))
    r_m = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample1 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"r"])
    r_m_mean = mean(r_m)
    auc_m = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample1 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"auc_l"])
    auc_m_mean = mean(auc_m)
    num_rep = length(r_m)
    # getting growth rates for control sample
    for(c in 1:length(controls)) {
      control = controls[c]
      sample2 = controls_strain_ids[c]
      print(paste(c, sample2, sep = " "))
      r_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"r"])
      r_w_mean = mean(r_w)
      auc_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"auc_l"])
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
      results_table = rbind(results_table, c(knockout, nitrogen_source, knockout, control, r_m_mean, r_w_mean, r_p_value_wt, r_p_value_aov, auc_m_mean, auc_w_mean, auc_p_value_wt, auc_p_value_aov, num_rep))
    }
}
}
results_table = results_table[-1,]
output1 = "growth_parameters_comparisons.transposons_vs_controls.csv";
write.table(results_table, file = output1, sep = "\t", row.names = F, col.names = T, quote = F)



# 2. comparing fitted growth parameters of nasD mutants against nasD wildtype

# controls (JE2 and spa) across all nitrogen sources

mutants = c("14324_8#22_nasD_p.Glu246Gln_mutant", "9716_3#67_nasD_p.Cys452Ser_mutant","8525_1#72_nasD_p.Thr656Ile_mutant")
mutants_strain_ids = c("146083", "58870","80131")
controls = c("14200_8#21_nasD_p.Glu246Gln_wildtype", "14324_8#72_nasD_p.Cys452Ser_wildtype","14208_8#20_nasD_p.Thr656Ile_wildtype")
controls_strain_ids = c("103356", "74975","144398")
mutations = c("nasD_p.Glu246Gln", "nasD_p.Cys452Ser", "nasD_p.Thr656Ile")

results_table = mat.or.vec(1,13)
# NOTE: the p-values of one-way ANOVA and wilcox tests to compare means are saved, both for growth rate comparison (r) and AUC comparison (AUC)
colnames(results_table) = c("mutation","nitrogen_source","sample1","sample2","r1_mean","r2_mean","r_wt_p_value","r_aov_p_value","auc1_mean","auc2_mean","auc_wt_p_value","auc_aov_p_value","num_rep")

#controls = c("nasD_knockout","spa_knockout", "JE2_control")
#controls_strain_ids = c("NE857","NE286", "JE2")

nitrogen_sources = unique(sort(metadata$nitrogen_source_pm3))

for(s in 1:length(nitrogen_sources))
{
  nitrogen_source = nitrogen_sources[s]
  print(paste(s, nitrogen_source, sep = " "))
  
  for(i in 1:length(mutants))
  {
    mutant = mutants[i]
    print(paste(i, mutant, sep = " "))
    
    # getting growth paramters for knockout sample
    sample1 = mutants_strain_ids[i]
    print(paste(i, sample1, sep = " "))
    r_m = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample1 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"r"])
    r_m_mean = mean(r_m)
    auc_m = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample1 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"auc_l"])
    auc_m_mean = mean(auc_m)
    num_rep = length(r_m)
    # getting growth rates for control sample
    for(c in 1:length(controls)) {
      if(c == i){
      control = controls[c]
      sample2 = controls_strain_ids[c]
      print(paste(c, sample2, sep = " "))
      r_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"r"])
      r_w_mean = mean(r_w)
      auc_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"auc_l"])
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
      results_table = rbind(results_table, c(mutations[i], nitrogen_source, mutant, control, r_m_mean, r_w_mean, r_p_value_wt, r_p_value_aov, auc_m_mean, auc_w_mean, auc_p_value_wt, auc_p_value_aov, num_rep))
      }
    }
  }
}

# 3. comparing fitted growth parameters of nasD mutants against controls


mutants = c("14324_8#22_nasD_p.Glu246Gln_mutant", "9716_3#67_nasD_p.Cys452Ser_mutant","8525_1#72_nasD_p.Thr656Ile_mutant")
mutants_strain_ids = c("146083", "58870","80131")
mutations = c("nasD_p.Glu246Gln", "nasD_p.Cys452Ser", "nasD_p.Thr656Ile")
controls = c("nasD_knockout","spa_knockout", "JE2_control")
controls_strain_ids = c("NE857","NE286", "JE2")

for(s in 1:length(nitrogen_sources))
{
  nitrogen_source = nitrogen_sources[s]
  print(paste(s, nitrogen_source, sep = " "))
  
  for(i in 1:length(mutants))
  {
    mutant = mutants[i]
    print(paste(i, mutant, sep = " "))
    
    # getting growth paramters for knockout sample
    sample1 = mutants_strain_ids[i]
    print(paste(i, sample1, sep = " "))
    r_m = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample1 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"r"])
    r_m_mean = mean(r_m)
    auc_m = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample1 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"auc_l"])
    auc_m_mean = mean(auc_m)
    num_rep = length(r_m)
    # getting growth rates for control sample
    for(c in 1:length(controls)) {
        control = controls[c]
        sample2 = controls_strain_ids[c]
        print(paste(c, sample2, sep = " "))
        r_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"r"])
        r_w_mean = mean(r_w)
        auc_w = as.numeric(curve_growth_par[which(curve_growth_par[,"strain_id"] == sample2 & curve_growth_par[,"nitrogen_source_pm3"] == nitrogen_source),"auc_l"])
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
        results_table = rbind(results_table, c(mutations[i], nitrogen_source, mutant, control, r_m_mean, r_w_mean, r_p_value_wt, r_p_value_aov, auc_m_mean, auc_w_mean, auc_p_value_wt, auc_p_value_aov, num_rep))
    }
  }
}

results_table = results_table[-1,]
output2 = "growth_parameters_comparisons.nasD_mutants_vs_controls.csv";
write.table(results_table, file = output2, sep = "\t", row.names = F, col.names = T, quote = F)





