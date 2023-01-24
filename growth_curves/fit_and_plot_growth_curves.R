# This R script is used to analyse microbial growth curves using growthcurver as done in:
# https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html
# https://rpubs.com/angelov/growthcurver
# And to plot fitted growth curved and extract growth parameters

## Required R libraries

library(dplyr)
library(reshape2)
library(ggplot2)
# install.packages("growthcurver")
library(growthcurver)
library(purrr)


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
metadata$label = paste(metadata$sample_id, metadata$gene, metadata$mutation, metadata$type, sep = "_")
metadata$mutation2 = paste(metadata$gene, metadata$mutation, sep = "_")


## Averaging values across replicates, only under daptomycin
unique_samples = unique(as.vector(metadata$label))
length(unique_samples)
# [1] 35

avgdata = mat.or.vec(nrow(rawdata), length(unique_samples) + 1)
colnames(avgdata) = c("time",unique_samples)
avgdata[,"time"] = rawdata[,"time"]

for(s in 1:length(unique_samples))
{
  sample = unique_samples[s]
  well_ids = metadata$curve_id[which(metadata$label==sample & metadata$environment==19)]
  avg = apply(rawdata[,well_ids], 1, mean)
  avgdata[,sample] = avg
}

avgdata = as.data.frame(avgdata)

## Averaging values across replicates, only no daptomycin

avgdata0 = mat.or.vec(nrow(rawdata), length(unique_samples) + 1)
colnames(avgdata0) = c("time",unique_samples)
avgdata0[,"time"] = rawdata[,"time"]

for(s in 1:length(unique_samples))
{
  sample = unique_samples[s]
  well_ids = metadata$curve_id[which(metadata$label==sample & metadata$environment==0)]
  avg = apply(rawdata[,well_ids], 1, mean)
  avgdata0[,sample] = avg
}
avgdata0 = as.data.frame(avgdata0)



## Fitting growth curves and plotting
## and extracting growth parameters

# ordering samples by gene and mutation
samples = colnames(avgdata0)[-1] # removing first "time" column
samples = samples[-length(samples)] # removing last NEG column
tmp = match(samples, metadata$label)
gene_mut_label = cbind(metadata$gene[tmp], metadata$mutation2[tmp], samples)
gene_mut_label = gene_mut_label[order(gene_mut_label[,1], gene_mut_label[,2]),]
samples = as.vector(gene_mut_label[,3])

# creating table to save growth parameters
growth_param0 = mat.or.vec(length(samples), 6)
colnames(growth_param0) = c("growth_rate","growth_rate_se","growth_rate_CI_left","growth_rate_CI_right","auc_l", "auc_e")

# NOTE: formulate to calculate CI extracted from here:
# https://stats.stackexchange.com/questions/474140/confidence-intervals-for-the-parameters-of-a-logistic-growth-curve

# no daptomycin present

# plotting 
pdf(file="fitted_growth_curves.no_dap.pdf")
par(mfcol = c(7,5)) # at least 35 slots are needed
par(mar=c(2, 2, 2, 2))
for(s in 1:length(samples))
{
  print(s)
  sample = samples[s]
  nas = which(is.na(avgdata0[,sample]))
  # making sure no NAs are available
  if(length(nas)==0) {
  # Fitting growth curve
  gc_fit <- SummarizeGrowth(avgdata0$time, avgdata0[,sample])
  # Saving growth parameters
  growth_rate = gc_fit$vals$r; growth_rate_se = gc_fit$vals$r_se;
  growth_rate_CI_left = growth_rate - qt(0.975, gc_fit$vals$df) * growth_rate_se;
  growth_rate_CI_right = growth_rate + qt(0.975, gc_fit$vals$df) * growth_rate_se;
  auc_l = gc_fit$vals$auc_l; auc_e = gc_fit$vals$auc_e;
  parameters_s = c(growth_rate, growth_rate_se, growth_rate_CI_left, growth_rate_CI_right, auc_l, auc_e)
  growth_param0[s,] = parameters_s
  
  # Plotting data
  title = sample; title = gsub("knockout_knockout", "knockout", title)
  plot(gc_fit$data$t, gc_fit$data$N, main = title, pch = 20, cex.main = 0.6, cex.axis=0.5, ylim = c(0,2))
  lines(gc_fit$data$t, predict(gc_fit$model), col = "red")
  text(2, 1.5, paste("r=",format(round(growth_rate, 2), nsmall = 2),sep = ""), cex = 0.5)
  text(10, 1.5, paste("AUC=",format(round(auc_e, 2), nsmall = 2),sep = ""), cex = 0.5)
  }
}
dev.off()

# saving growth parameters
colnames(gene_mut_label) = c("gene", "mutation", "sample")
out_table = cbind(gene_mut_label, growth_param0)
write.table(out_table, file = "fitted_growth_curves.no_dap.growth_parameters.csv", sep = "\t", col.names = T, row.names = F, quote = F)




# daptomycin present

# creating table to save growth parameters
growth_param = mat.or.vec(length(samples), 6)
colnames(growth_param) = c("growth_rate","growth_rate_se","growth_rate_CI_left","growth_rate_CI_right","auc_l", "auc_e")

pdf(file="fitted_growth_curves.dap.pdf")
par(mfcol = c(7,5)) # at least 35 slots are needed
par(mar=c(2, 2, 2, 2))
for(s in 1:length(samples))
{
  print(s)
  sample = samples[s]
  nas = which(is.na(avgdata[,sample]))
  # making sure no NAs are available
  if(length(nas)==0) {
    # fitting the curve
    gc_fit <- SummarizeGrowth(avgdata$time, avgdata[,sample])
    # Saving growth parameters
    growth_rate = gc_fit$vals$r; growth_rate_se = gc_fit$vals$r_se;
    growth_rate_CI_left = growth_rate - qt(0.975, gc_fit$vals$df) * growth_rate_se;
    growth_rate_CI_right = growth_rate + qt(0.975, gc_fit$vals$df) * growth_rate_se;
    auc_l = gc_fit$vals$auc_l; auc_e = gc_fit$vals$auc_e;
    parameters_s = c(growth_rate, growth_rate_se, growth_rate_CI_left, growth_rate_CI_right, auc_l, auc_e)
    growth_param[s,] = parameters_s
    # editing plot title
    title = sample; title = gsub("knockout_knockout", "knockout", title)
    plot(gc_fit$data$t, gc_fit$data$N, main = title, pch = 20, cex.main = 0.6, cex.axis=0.5, ylim = c(0,2))
    if(gc_fit$vals$k>0) {
      lines(gc_fit$data$t, predict(gc_fit$model), col = "red")
      text(2, 1.5, paste("r=",format(round(growth_rate, 2), nsmall = 2),sep = ""), cex = 0.5)
      text(10, 1.5, paste("AUC=",format(round(auc_e, 2), nsmall = 2),sep = ""), cex = 0.5)
    }
  }
}
dev.off()

# saving growth parameters
colnames(gene_mut_label) = c("gene", "mutation", "sample")
out_table = cbind(gene_mut_label, growth_param)
write.table(out_table, file = "fitted_growth_curves.dap.growth_parameters.csv", sep = "\t", col.names = T, row.names = F, quote = F)




# daptomycin present - fitting after 7 hours

# creating table to save growth parameters
growth_param7h = mat.or.vec(length(samples), 6)
colnames(growth_param7h) = c("growth_rate","growth_rate_se","growth_rate_CI_left","growth_rate_CI_right","auc_l", "auc_e")

pdf(file="fitted_growth_curves.dap.after7h.pdf")
par(mfcol = c(7,5)) # at least 35 slots are needed
par(mar=c(2, 2, 2, 2))
for(s in 1:length(samples))
{
  print(s)
  sample = samples[s]
  nas = which(is.na(avgdata[,sample]))
  # making sure no NAs are available
  if(length(nas)==0) {
    # Fitting the curve
    subset = which(avgdata$time >= 7)
    gc_fit0 <- SummarizeGrowth(avgdata$time, avgdata[,sample])
    gc_fit7 <- SummarizeGrowth(avgdata$time[subset], avgdata[subset,sample])
    # Saving growth parameters
    growth_rate = gc_fit7$vals$r; growth_rate_se = gc_fit7$vals$r_se;
    growth_rate_CI_left = growth_rate - qt(0.975, gc_fit7$vals$df) * growth_rate_se;
    growth_rate_CI_right = growth_rate + qt(0.975, gc_fit7$vals$df) * growth_rate_se;
    auc_l = gc_fit7$vals$auc_l; auc_e = gc_fit7$vals$auc_e;
    parameters_s = c(growth_rate, growth_rate_se, growth_rate_CI_left, growth_rate_CI_right, auc_l, auc_e)
    growth_param7h[s,] = parameters_s
    # plotting curves
    title = sample; title = gsub("knockout_knockout", "knockout", title)
    plot(gc_fit0$data$t, gc_fit0$data$N, main = title, pch = 20, cex.main = 0.6, cex.axis=0.5, ylim = c(0,2))
    if(gc_fit7$vals$k>0) {
      lines(gc_fit7$data$t, predict(gc_fit7$model), col = "red")
      text(2, 1.5, paste("r=",format(round(growth_rate, 2), nsmall = 2),sep = ""), cex = 0.5)
      text(10, 1.5, paste("AUC=",format(round(auc_e, 2), nsmall = 2),sep = ""), cex = 0.5)
    }
  }
}
dev.off()

# saving growth parameters
colnames(gene_mut_label) = c("gene", "mutation", "sample")
out_table = cbind(gene_mut_label, growth_param7h)
write.table(out_table, file = "fitted_growth_curves.dap.after7h.growth_parameters.csv", sep = "\t", col.names = T, row.names = F, quote = F)




# daptomycin present - fitting after 10 hours

# creating table to save growth parameters
growth_param10h = mat.or.vec(length(samples), 6)
colnames(growth_param10h) = c("growth_rate","growth_rate_se","growth_rate_CI_left","growth_rate_CI_right","auc_l", "auc_e")

pdf(file="fitted_growth_curves.dap.after10h.pdf")
par(mfcol = c(7,5)) # at least 35 slots are needed
par(mar=c(2, 2, 2, 2))
for(s in 1:length(samples))
{
  print(s)
  sample = samples[s]
  nas = which(is.na(avgdata[,sample]))
  # making sure no NAs are available
  if(length(nas)==0) {
    # Fitting the curve
    subset = which(avgdata$time >= 10)
    gc_fit0 <- SummarizeGrowth(avgdata$time, avgdata[,sample])
    gc_fit10 <- SummarizeGrowth(avgdata$time[subset], avgdata[subset,sample])
    # Saving growth parameters
    growth_rate = gc_fit10$vals$r; growth_rate_se = gc_fit10$vals$r_se;
    growth_rate_CI_left = growth_rate - qt(0.975, gc_fit10$vals$df) * growth_rate_se;
    growth_rate_CI_right = growth_rate + qt(0.975, gc_fit10$vals$df) * growth_rate_se;
    auc_l = gc_fit10$vals$auc_l; auc_e = gc_fit10$vals$auc_e;
    parameters_s = c(growth_rate, growth_rate_se, growth_rate_CI_left, growth_rate_CI_right, auc_l, auc_e)
    growth_param7h[s,] = parameters_s
    # plotting curves
    title = sample; title = gsub("knockout_knockout", "knockout", title)
    plot(gc_fit0$data$t, gc_fit0$data$N, main = title, pch = 20, cex.main = 0.6, cex.axis=0.5, ylim = c(0,2))
    if(gc_fit10$vals$k>0) {
      lines(gc_fit10$data$t, predict(gc_fit10$model), col = "red")
      text(2, 1.5, paste("r=",format(round(growth_rate, 2), nsmall = 2),sep = ""), cex = 0.5)
      text(10, 1.5, paste("AUC=",format(round(auc_e, 2), nsmall = 2),sep = ""), cex = 0.5)
    }
  }
}
dev.off()

# saving growth parameters
colnames(gene_mut_label) = c("gene", "mutation", "sample")
out_table = cbind(gene_mut_label, growth_param10h)
write.table(out_table, file = "fitted_growth_curves.dap.after10h.growth_parameters.csv", sep = "\t", col.names = T, row.names = F, quote = F)


