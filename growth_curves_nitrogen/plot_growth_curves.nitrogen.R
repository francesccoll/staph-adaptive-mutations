
# This R script is used to plot microbial growth curves as done in:
# http://bconnelly.net/posts/analyzing_microbial_growth_with_r/

library(reshape2)
library(dplyr)
library(ggplot2)

## Reading and merging input data

rawdata = read.delim("Export_All_Plates_Single_File.curves_growthdata.csv", sep = "\t", header = T)
dim(rawdata)
# [1]  193 2881

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

# Reshape the data. Instead of rows containing the time and OD600 readings for each Well, 
# rows will contain the time, a Well ID (curve_id), and the reading at that Well.
reshaped <- melt(rawdata, id=c("time"), variable.name="curve_id", value.name="OD600")
dim(reshaped)
# [1] 555840      3
# Fixing curve_id
reshaped$curve_id = gsub("^X", "", reshaped$curve_id)

# Add information about the strain from the metadata table. For each curve_id/well
# defined in both the reshaped data and metadata, each resulting row
# will contain the absorbance measurement as well as the additional columns
# and values from the metadata
annotated <- inner_join(reshaped, metadata, by="curve_id")
dim(annotated)
# [1] 555840     15

# Save the annotated data as a CSV for storing, sharing, etc.
write.table(annotated, file = "Export_All_Plates_Single_File.annotated.csv", sep = "\t", col.names = T, row.names = F, quote = F)

conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}


# Group the data by the different experimental variables and calculate the
# sample size, average OD600, and 95% confidence limits around the mean
# among the replicates. Also remove all records where the Strain is NA.
# stats <- annotated %>%
#   group_by(nitrogen_source_pm3, label, time, gene, mutation2) %>%
#   summarise(N=length(OD600),
#             Average=mean(OD600),
#             CI95=conf_int95(OD600)) %>%
#   filter(!is.na(label))

stats <- annotated %>%
  group_by(nitrogen_source_pm3, label, time, gene, mutation2) %>%
  summarise(N=length(OD600),
            Average=mean(OD600),
            sd=sd(OD600)) %>%
  filter(!is.na(label))


dim(stats)
# [1] 185280      8

head(stats)
# # A tibble: 6 × 8
# # Groups:   nitrogen_source_pm3, label, time, gene [6]
# nitrogen_source_pm3 label                                 time gene  mutation2            N Average  CI95
# <chr>               <chr>                                <dbl> <chr> <chr>            <int>   <dbl> <dbl>
# 1 Acetamide           14200_8#21_nasD_p.Glu246Gln_wildtype  0    nasD  nasD_p.Glu246Gln     3    21.4  5.25
# 2 Acetamide           14200_8#21_nasD_p.Glu246Gln_wildtype  0.25 nasD  nasD_p.Glu246Gln     3    23.1 10.3 
# 3 Acetamide           14200_8#21_nasD_p.Glu246Gln_wildtype  0.5  nasD  nasD_p.Glu246Gln     3    19.9  3.33
# 4 Acetamide           14200_8#21_nasD_p.Glu246Gln_wildtype  0.75 nasD  nasD_p.Glu246Gln     3    19.3  4.40
# 5 Acetamide           14200_8#21_nasD_p.Glu246Gln_wildtype  1    nasD  nasD_p.Glu246Gln     3    22.6  5.55
# 6 Acetamide           14200_8#21_nasD_p.Glu246Gln_wildtype  1.25 nasD  nasD_p.Glu246Gln     3    21.6  4.61



# nitrogen gene knockouts + wildtype strain

nitrogen_genes = c("nasD", "ureG", "spa")

stats_mut = filter(stats, gene %in% nitrogen_genes & grepl("knockout",label))
stats_mut_control = filter(stats, grepl("JE2", label))
stats_mut = rbind(stats_mut, stats_mut_control)
stats_mut$label = gsub("_knockout_knockout", "_knockout", stats_mut$label)
stats_mut$label = gsub("_none_none", "_none", stats_mut$label)

nitrogen_sources = c("Negative_Control", "Urea", "Nitrate", "Nitrite", "Ammonia", "L-Glutamine", "L-Glutamic_Acid","Glycine","L-Methionine","D-Valine","L-Glutamine")
nitrogen_sources = unique(metadata$nitrogen_source_pm3)

for(s in 1:length(nitrogen_sources))
{
  nitrogen_source = nitrogen_sources[s]
  print(nitrogen_source)
  plot_title = paste("nitrogen gene knockouts", nitrogen_source, sep=" ")
  
  stats_env = filter(stats_mut, nitrogen_source_pm3 == nitrogen_source)
  
  gp = ggplot(data=stats_env, aes(x=time, y=Average, color=label)) +
    geom_ribbon(aes(ymin=Average-sd, ymax=Average+sd, fill=label), color=NA, alpha=0.3) +
    # geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=label), color=NA, alpha=0.3) +
    geom_line() +
    # scale_y_log10() +
    facet_grid(nitrogen_source_pm3 ~ .) + ylim(1, 250) +
    # theme(legend.position="bottom", legend.box = "vertical") +
    labs(x="Time (Hours)", y="Absorbance at 600 nm", title = plot_title)

plot_file = paste("nitrogen_genes_knockouts.",nitrogen_source,".growth_curves.pdf",sep="")
ggsave(plot_file, plot = gp, device = "pdf", width = 8, height = 4, dpi = 300, units = "in")
}


# nasD mutations + controls
mutations = c("nasD_p.Glu246Gln", "nasD_p.Thr656Ile", "nasD_p.Cys452Ser")
nitrogen_sources = c("Negative_Control", "Urea", "Nitrate", "Nitrite", "Ammonia", "L-Glutamine", "L-Glutamic_Acid","Glycine","L-Methionine","D-Valine","L-Glutamine","L-Leucine","L-Isoleucine","L-Tryptophan")
nitrogen_sources = unique(metadata$nitrogen_source_pm3)

for(s in 1:length(nitrogen_sources))
{
  nitrogen_source = nitrogen_sources[s]
  print(nitrogen_source)
  
for(m in 1:length(mutations))
{
  # choosing controls to include
  gene = unlist(strsplit(mutations[m], "_"))[1]
  knockout = paste(gene, "_knockout", sep = "")
  print(mutations[m])
  
  # choosing strains to plot
  # stats_mut = filter(stats, mutation2 == mutations[m] | mutation2 == knockout)
  # adding transposon wildtype control strain too
  stats_env = filter(stats, nitrogen_source_pm3 == nitrogen_source)
  stats_mut = filter(stats_env, mutation2 == mutations[m] | mutation2 == knockout | grepl("JE2", label))
  plot_title = paste(mutations[m], nitrogen_source, sep=" ")
  gp = ggplot(data=stats_mut, aes(x=time, y=Average, color=label)) +
    # geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=label), color=NA, alpha=0.3) + ylim(1, 250) +
    geom_ribbon(aes(ymin=Average-sd, ymax=Average+sd, fill=label), color=NA, alpha=0.3) + ylim(1, 250) +
    geom_line() + facet_grid(nitrogen_source_pm3 ~ .) + labs(x="Time (Hours)", y="Absorbance at 600 nm", title = plot_title)
  
  plot_file = paste(mutations[m], ".", nitrogen_source, ".growth_curves.pdf", sep = "")
  ggsave(plot_file, plot = gp, device = "pdf", width = 8, height = 4, dpi = 300, units = "in")
}
}


### Plotting single growth curve replicates for sample nasD_p.Thr656Ile (the wildtype had wider CI and SD)


# nasD mutations + controls
mutations = c("nasD_p.Thr656Ile")
nitrogen_sources = c("Negative_Control", "Urea", "Nitrite", "Ammonia")
  
for(s in 1:length(nitrogen_sources))
{
  nitrogen_source = nitrogen_sources[s]
  print(nitrogen_source)
  
  for(m in 1:length(mutations))
  {
    # choosing controls to include
    gene = unlist(strsplit(mutations[m], "_"))[1]
    knockout = paste(gene, "_knockout", sep = "")
    print(mutations[m])
    
    # choosing strains to plot
    # adding transposon wildtype control strain too
    annotated_env = filter(annotated, nitrogen_source_pm3 == nitrogen_source)
    annotated_mut = filter(annotated_env, mutation2 == mutations[m])
    plot_title = paste(mutations[m], "replicates",  nitrogen_source,sep=" ")
    gp = ggplot(data=annotated_mut, aes(x=time, y=OD600, group = interaction(type, repeat_id), color=repeat_id, linetype=type)) +
      # geom_ribbon(aes(ymin=Average-sd, ymax=Average+sd, fill=label), color=NA, alpha=0.3) + 
      ylim(1, 250) +
      geom_line() + facet_grid(nitrogen_source_pm3 ~ .) + labs(x="Time (Hours)", y="Absorbance at 600 nm", title = plot_title)
    
    plot_file = paste(mutations[m], ".", nitrogen_source, ".growth_curves.replicates.pdf", sep = "")
    ggsave(plot_file, plot = gp, device = "pdf", width = 8, height = 4, dpi = 300, units = "in")
  }
}



