# This R script is used to summarise within-host genetic diversity across all patients

require(gdata)
require(ggplot2)

#########################################################################################################################################################
#####                                                   Creating summary table of per patient diversity                                           #######
#########################################################################################################################################################


metadata_file = "/project/directory/data/colonising_pairs_patient-matched_gwas_dataset.all_pairs.metadata.xlsx";
lmT = read.xls(metadata_file, sheet = 1, header = T)
dim(lmT)
# [1] 2627   22

mutations_file = "/project/directory/data/colonising_pairs_patient-matched_gwas_dataset.mrca.snippy_snpEff_annotated_variants.nctc8325.in_ref_VCF.filtered.txt.gz";
lmM = read.csv(gzfile(mutations_file), sep = "\t", header = T)
dim(lmM)
# [1] 5552   32

# keep only QCed pairs
lmT = subset(lmT, QC1 == "yes")
dim(lmT)
# [1] 2270   22

collections = unique(as.vector(lmT$collection))

lmSummary = mat.or.vec(1,10)
colnames(lmSummary) = c("collection", "patient", "num_pairs", "isolates_num", "isolates_text", "min_core_snps", "max_core_snps", "min_time_span", "max_time_span", "num_wVar")

for(c in 1:length(collections))
{
  collection = as.character(collections[c])
  patients = unique(as.vector(lmT$patientID[which(lmT$collection==collection)]))
  for(p in 1:length(patients))
  {
    patient = as.character(patients[p])
    tmp = which(lmT$patientID == patient)
    num_pairs = as.character(length(tmp))
    isolates = unique(c(as.vector(lmT$isolateID1[tmp]), as.vector(lmT$isolateID2[tmp])))
    isolates_text = paste(isolates, collapse = ";")
    isolates_num = as.character(length(isolates))
    core_snps = as.vector(lmT$SNPdisCoreGenome[tmp])
    min_core_snps = as.character(min(core_snps)); max_core_snps = as.character(max(core_snps));
    # Time distance between isolates
    max_time_span = "not_calculated"; min_time_span = "not_calculated";
    if(collection=="paterson2015" | collection=="chow2017"){ max_time_span = "unknown"; min_time_span = "unknown"; }
    else{ time_distances = as.Date(lmT$collection_date2[tmp], format = "%Y-%m-%d") - as.Date(lmT$collection_date1[tmp], format = "%Y-%m-%d"); max_time_span = as.character(max(as.integer(time_distances))); min_time_span = as.character(min(as.integer(time_distances))); }
    # Extracting number of whole-genome mutations
    tmp2 = which(lmM$patient_id == patient)
    num_wVar = as.character(length(tmp2))
    # Saving information
    newrow = c(collection, patient, num_pairs, isolates_num, isolates_text, min_core_snps, max_core_snps, min_time_span, max_time_span, num_wVar)
    lmSummary = rbind(lmSummary, newrow)
  }
}
lmSummary = lmSummary[-1,]
dim(lmSummary)
# [1] 791  10

# 791 number of patients after filtering

sum(as.integer(lmSummary[,"isolates_num"]))
# [1] 3060

summary_table_file = "colonising_pairs_patient-matched_gwas_dataset.all_pairs.diversity.csv";
write.table(lmSummary, file=summary_table_file, sep = '\t', col.names = T, row.names = F, quote = F)



#########################################################################################################################################################
#####                                                       Creating per patient diversity barplots                                              #######
#########################################################################################################################################################


lmSummary = read.delim(summary_table_file, sep = "\t", header = T)
dim(lmSummary)
# [1] 791  10

size_dot = 1; size_axis_lines = 0.3; text_y_offset = 4; font = "Times"; dot_color = "dimgray";
axis_text_size = 15; axis_title_size = 15; ann_text_size = 5;
g1 <- ggplot(lmSummary, aes(x=max_core_snps)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100), width = 0.8) + 
  xlim(-1, 51) + 
  ylim(0, 40) + 
  xlab("Number of core-genome SNPs") +
  ylab("Percentage of Individuals") +
  ggtitle("Within-host genetic diversity") +
  theme(text = element_text(family = font)) +
  theme(axis.text = element_text(size=axis_text_size, color="black"), axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))
plot_file = "within-host_diversity.core-genome_SNPs.pdf";
ggsave(plot_file, plot = g1, device = "pdf",width = 8, height = 5, dpi = 300, units = "in")




# If num_wVar > 50, then make it 50 so that the last bar represents individuals with >=50 SNPs
lmSummary[which(lmSummary$num_wVar > 50),"num_wVar"] = 50

size_dot = 1; size_axis_lines = 0.3; text_y_offset = 4; font = "Times"; dot_color = "dimgray";
axis_text_size = 15; axis_title_size = 15; ann_text_size = 5;
g1 <- ggplot(lmSummary, aes(x=num_wVar)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100), width = 0.8) + 
  xlim(-1, 51) + 
  ylim(0, 40) + 
  xlab("Number of whole-genome variants") +
  ylab("Percentage of Individuals") +
  ggtitle("Within-host genetic diversity") +
  theme(text = element_text(family = font)) +
  theme(axis.text = element_text(size=axis_text_size, color="black"), axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))
plot_file = "within-host_diversity.whole-genome_variants.pdf";
ggsave(plot_file, plot = g1, device = "pdf",width = 8, height = 5, dpi = 300, units = "in")


















