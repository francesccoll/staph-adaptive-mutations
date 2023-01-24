
# This R script is used to plot microbial growth curves as done in:
# http://bconnelly.net/posts/analyzing_microbial_growth_with_r/

library(reshape2)
library(dplyr)
library(ggplot2)

## Reading and merging input data

rawdata = read.delim("curves_all_14-10-22.csv", sep = ",", header = T)
dim(rawdata)
# [1]  48 598
colnames(rawdata)[1] = "time"

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

# Reshape the data. Instead of rows containing the time and OD600 readings for each Well, 
# rows will contain the time, a Well ID (curve_id), and the reading at that Well.
reshaped <- melt(rawdata, id=c("time"), variable.name="curve_id", value.name="OD600")
dim(reshaped)
# [1] 28656     3


# Add information about the strain from the metadata table. For each curve_id/well
# defined in both the reshaped data and metadata, each resulting row
# will contain the absorbance measurement as well as the additional columns
# and values from the metadata
annotated <- inner_join(reshaped, metadata, by="curve_id")
dim(annotated)
# [1] 28656    11

# Save the annotated data as a CSV for storing, sharing, etc.
write.table(annotated, file = "curves_all_14-10-22.annotated.csv", sep = "\t", col.names = T, row.names = F, quote = F)

conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}


# Group the data by the different experimental variables and calculate the
# sample size, average OD600, and 95% confidence limits around the mean
# among the replicates. Also remove all records where the Strain is NA.
stats <- annotated %>%
  group_by(environment, label, time, gene, mutation2) %>%
  summarise(N=length(OD600),
            Average=mean(OD600),
            CI95=conf_int95(OD600)) %>%
  filter(!is.na(label))

dim(stats)
# [1] 3216    8

head(stats)
# # A tibble: 6 × 8
# # Groups:   environment, label, time, gene [6]
# environment label                              time gene  mutation2           N Average   CI95
# <chr>       <chr>                             <dbl> <chr> <chr>           <int>   <dbl>  <dbl>
# 1 0           10770_3#44_vraA_p.Lys419fs_mutant   0   vraA  vraA_p.Lys419fs     9   0.260 0.0671
# 2 0           10770_3#44_vraA_p.Lys419fs_mutant   0.5 vraA  vraA_p.Lys419fs     9   0.235 0.0516
# 3 0           10770_3#44_vraA_p.Lys419fs_mutant   1   vraA  vraA_p.Lys419fs     9   0.220 0.0465
# 4 0           10770_3#44_vraA_p.Lys419fs_mutant   1.5 vraA  vraA_p.Lys419fs     9   0.190 0.0406
# 5 0           10770_3#44_vraA_p.Lys419fs_mutant   2   vraA  vraA_p.Lys419fs     9   0.215 0.0436
# 6 0           10770_3#44_vraA_p.Lys419fs_mutant   2.5 vraA  vraA_p.Lys419fs     9   0.211 0.0434


## Plot the average OD600 over time for each strain in each environment

# First plot, virulence factors + wildtype strain
virulence_factors = c("clfA", "clfB", "fnbA", "fnbB", "spa")
stats_mut = filter(stats, gene %in% virulence_factors | grepl("JE2",label))

gp = ggplot(data=stats_mut, aes(x=time, y=Average, color=label)) +
  geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=label),
              color=NA, alpha=0.3) + 
  geom_line() +
  # scale_y_log10() +
  facet_grid(environment ~ .) +
  # theme(legend.position="bottom", legend.box = "vertical") + 
  labs(x="Time (Hours)", y="Absorbance at 600 nm")

plot_file = "virulence_factors.knockouts.growth_curves.pdf"
ggsave(plot_file, plot = gp, device = "pdf", width = 8, height = 4, dpi = 300, units = "in")


# Second plot, daptomycin resistance genes knockouts + wildtype strain
daptomycin_genes = c("mprF", "clpP", "pitA", "vraA", "vraB", "proP", "pstS","pstC","cvfB")
daptomycin_genes = c("mprF", "clpP", "pitA", "vraA", "pstS")

stats_mut = filter(stats, gene %in% daptomycin_genes & grepl("knockout",label))
stats_mut_control = filter(stats, grepl("JE2", label))
stats_mut = rbind(stats_mut, stats_mut_control)
  
gp = ggplot(data=stats_mut, aes(x=time, y=Average, color=label)) +
  geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=label),
              color=NA, alpha=0.3) + 
  geom_line() +
  # scale_y_log10() +
  facet_grid(environment ~ .) + 
  # theme(legend.position="bottom", legend.box = "vertical") + 
  labs(x="Time (Hours)", y="Absorbance at 600 nm")

plot_file = "daptomycin_genes.knockouts.growth_curves.pdf"
ggsave(plot_file, plot = gp, device = "pdf", width = 8, height = 4, dpi = 300, units = "in")


# vraA and pstS mutations + controls
mutations = c("vraA_p.Lys419fs", "vraA_p.Pro60Ser", "vraA_p.Gln13", "vraA_p.Leu417","pstS_p.Ser164fs","pstS_p.Gly206Glu","pstS_p.Gln217")

for(m in 1:length(mutations))
{
  # choosing controls to include
  gene = unlist(strsplit(mutations[m], "_"))[1]
  knockout = paste(gene, "_knockout", sep = "")

  # choosing strains to plot
  # stats_mut = filter(stats, mutation2 == mutations[m] | mutation2 == knockout)
  # adding transposon wildtype control strain too
  stats_mut = filter(stats, mutation2 == mutations[m] | mutation2 == knockout | grepl("JE2", label))
  gp = ggplot(data=stats_mut, aes(x=time, y=Average, color=label)) +
    geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=label),
                color=NA, alpha=0.3) + 
    geom_line() + facet_grid(environment ~ .) + labs(x="Time (Hours)", y="Absorbance at 600 nm")
  
  plot_file = paste(mutations[m], ".growth_curves.pdf", sep = "")
  ggsave(plot_file, plot = gp, device = "pdf", width = 8, height = 4, dpi = 300, units = "in")
}




