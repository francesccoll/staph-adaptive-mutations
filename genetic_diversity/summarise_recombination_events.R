
# This R script is used to study recombination events identified in script poisson_test.v4.NCTC8325.R

###############################################################################################################
####                                         INPUT/OUTPUT FILES                                            ####
###############################################################################################################

dir = "/project/directory/data/"

# Filtered mutations after excluding recombination, and repeated mutations from multiple pairs of the same patient
filtered_mutations_table = paste(dir, "colonising_pairs_patient-matched_gwas_dataset.mrca.snippy_snpEff_annotated_variants.nctc8325.in_ref_VCF.filtered.txt.gz", sep = "")

# Filtered recombination mutations, without removing repeated mutations from multiple pairs of the same patient
recombination_mutations_table = paste(dir, "colonising_pairs_patient-matched_gwas_dataset.mrca.snippy_snpEff_annotated_variants.nctc8325.in_ref_VCF.recombination.txt.gz", sep = "")

output = "colonising_pairs_patient-matched_gwas_dataset.recombination_events.csv"

lmMF = read.delim(gzfile(filtered_mutations_table), sep = "\t", header = T)
dim(lmMF)
# [1] 5552   32

lmMR = read.delim(gzfile(recombination_mutations_table), sep = "\t", header = T)
dim(lmMR)
# [1] 3922   32


###############################################################################################################
####             Removing repeated mutations from multiple pairs of the same patient                       ####
###############################################################################################################

lmMR = lmMR[order(lmMR$sample),]

patients = unique(as.vector(lmMR$patient_id))
length(patients)
# [1] 117

mutations_to_remove = vector();
for(p in 1:length(patients))
{
  tmp = which(lmMR$patient_id == patients[p])
  tmp2 = which(duplicated(lmMR$pos[tmp]))
  if(length(tmp2)>0){ mutations_to_remove = c(mutations_to_remove, tmp[tmp2]); }
}
length(mutations_to_remove)
#[1] 2201
lmMR = lmMR[-mutations_to_remove,];
dim(lmMR)
# [1] 1721   32

# NOTE: a total of 1721 variants are attributable to recombination across 117 hosts


###############################################################################################################
####                               Recombination variants per host                                      ####
###############################################################################################################


sort(table(lmMR$patient_id))

###############################################################################################################
####                               Assigning recombination event Ids                                      ####
###############################################################################################################

lmMR = lmMR[order(as.vector(lmMR$patient_id), as.numeric(lmMR$pos)),]
lmMR$recombination_event = "NA"
laPatientsRecomb = unique(as.vector(lmMR$patient_id))
length(laPatientsRecomb)
# [1] 117

window = 1000; # a window of 100 failed to remove recombination in phage region
for(s in 1:length(laPatientsRecomb))
{
  sss = which(lmMR$patient_id==as.character(laPatientsRecomb[s]))
  mutationsSam = lmMR[sss,];
  mutationsSam = mutationsSam[order(as.numeric(mutationsSam$pos)),]
  recombination_event = 1
  recom_idx_prev = vector()
  for(m in 1:length(sss))
  {
    print(m)
    if(mutationsSam$recombination_event[m] == "NA")
    {
      recom_left_idx = which((mutationsSam$pos < mutationsSam$pos[m]) & (mutationsSam$pos >= (mutationsSam$pos[m]-window))); # check if neighbouring mutations downstream
      recom_right_idx = which((mutationsSam$pos > mutationsSam$pos[m]) & (mutationsSam$pos <= (mutationsSam$pos[m]+window))); # check if neighbouring mutations upstream
      recom_idx = c(recom_left_idx, recom_right_idx)
      if(length(recom_idx)>0)
      {
        recom_idx = c(m, recom_idx)
        if(length(recom_idx_prev)==0)
        {
          lmMR$recombination_event[sss[recom_idx]] = as.character(recombination_event)
          mutationsSam$recombination_event[recom_idx] = as.character(recombination_event)
        } else
        {
          if(length(intersect(recom_idx, recom_idx_prev))>0)
          {
            # if current and previous recombination positions overlap, then assign the same recombination_event
            lmMR$recombination_event[sss[recom_idx]] = as.character(recombination_event)
            mutationsSam$recombination_event[recom_idx] = as.character(recombination_event)
          } else
          {
            # if not, then assign a new recombination_event
            recombination_event = recombination_event + 1
            lmMR$recombination_event[sss[recom_idx]] = as.character(recombination_event)
            mutationsSam$recombination_event[recom_idx] = as.character(recombination_event)
          }
        }
      } else
      {
        recombination_event = recombination_event + 1
      }
      recom_idx_prev = recom_idx
    }
  }
}


lmMR$recombination_event = paste(lmMR$patient_id, "-re", lmMR$recombination_event, sep="")
dim(lmMR)
# [1] 1721   33


###############################################################################################################
####                            Recombination within and outside hotspots                                 ####
###############################################################################################################

# coordinates of phages phi11, phi12 and phi13 in the NCTC8325 reference genomes were derived from blasting
# their FASTA sequences (AF424781.1, AF424782.1 and AF424783.1) against that of NCTC8325 

hotspots_start = c(1923408, 1463618, 2031911)
hotspots_end = c(1967011, 1508580, 2074632)
hotsplot_pos = vector()

recomb_hot = vector()  # to store recombination variant positions within hotspots
point_hot = vector()  # to store point mutation variant positions within hotspots
for(h in 1:length(hotspots_start))
{
  tmp = which(lmMR$pos >= hotspots_start[h] & lmMR$pos <= hotspots_end[h])
  recomb_hot = c(recomb_hot, tmp)
  tmp2 = which(lmMF$pos >= hotspots_start[h] & lmMF$pos <= hotspots_end[h])
  point_hot = c(point_hot, tmp2)
  hotsplot_pos = c(hotsplot_pos, seq(hotspots_start[h], hotspots_end[h], 1))
}
length(unique(recomb_hot))
# [1] 1367
length(unique(point_hot))
# [1] 0  > this is expected because mutations in phage regions were removed
length(unique(hotsplot_pos))
# [1] 131289







###############################################################################################################
####                          Summarising recombination events per patient                                 ####
###############################################################################################################


lmRecomEvents = mat.or.vec(1,8)
colnames(lmRecomEvents) = c("patient_id", "recombination_event", "start", "end", "length", "num_var", "loci", "samples")

for(p in 1:length(laPatientsRecomb))
{
  patient_id = as.character(laPatientsRecomb[p])
  recombination_events = unique(lmMR$recombination_event[which(lmMR$patient_id==patient_id)])
  for(re in 1:length(recombination_events))
  {
    re_idx = which(lmMR$recombination_event == recombination_events[re])
    start = as.character(lmMR$pos[re_idx[1]])
    end = as.character(lmMR$pos[re_idx[length(re_idx)]])
    length = as.character(as.numeric(end) - as.numeric(start) + 1)
    num_var = as.character(length(re_idx))
    loci = as.character(paste(unique(lmMR$locus_tag[re_idx]),collapse = ";"))
    samples = as.character(paste(unique(lmMR$sample[re_idx]),collapse = ";"))
    newrow = c(patient_id, recombination_events[re], start, end, length, num_var, loci, samples)
    lmRecomEvents = rbind(lmRecomEvents, newrow)
  }
}

lmRecomEvents = lmRecomEvents[-1,]

write.table(lmRecomEvents, file=output, sep = '\t', col.names = T, row.names = F, quote = F)


###############################################################################################################
####                                         Summary statistics                                              ####
###############################################################################################################

length(names(table(lmRecomEvents[,"patient_id"])))
# [1] 117

quantile(table(lmRecomEvents[,"patient_id"]))
# 0%  25%  50%  75% 100% 
# 1    1    1    3   25 

quantile(as.numeric(lmRecomEvents[,"length"]))
# 0%    25%    50%    75%   100% 
# 3.0   37.0  214.5  665.0 6334.0 

quantile(as.numeric(lmRecomEvents[,"num_var"]))
# 0%  25%  50%  75% 100% 
# 2    2    2    5   65 


###############################################################################################################
####                    Calculating # genomes and density of mutations with recombinations                 ####
###############################################################################################################

window_size = 2000
window_step = 1000
chr_length = 2821360

window_start = seq(1,chr_length,window_step)
window_end = window_start + window_size -1

window_recom_table = mat.or.vec(1,6)
colnames(window_recom_table) =  c("window_i", "window_start", "window_end", "window_size", "num_patients", "mean_density")

for(w in 1:length(window_start))
{
  pos = which(lmMR$pos >= window_start[w] & lmMR$pos < window_end[w])
  print(paste(w, " - ", window_start[w], " - ", window_end[w], " - ", length(pos), sep=""))
  num_patients = 0
  mean_density = 0
  if(length(pos) > 0)
  {
    patients = unique(as.vector(lmMR$patient_id[pos]))
    num_patients = length(patients)
    num_variant = vector(); # vector to store density of variants per patient in window w
    for(p in 1:length(patients))
    {
      pos_pt = which(lmMR$pos >= window_start[w] & lmMR$pos < window_end[w] & lmMR$patient_id == patients[p])
      num_variant = c(num_variant, length(pos_pt))
    }
    density_variants_per_kp = (num_variant/window_size)*1000
    mean_density = mean(density_variants_per_kp)
  }
  newrow = c(as.character(w), as.character(window_start[w]), as.character(window_end[w]), as.character(window_size), as.character(num_patients), as.character(mean_density))
  window_recom_table = rbind(window_recom_table, newrow)
}
window_recom_table = window_recom_table[-1,]
dim(window_recom_table)
# [1] 2822    6

output2 = "colonising_pairs_patient-matched_gwas_dataset.recombination_density.csv"
write.table(window_recom_table, file=output2, sep = '\t', col.names = T, row.names = F, quote = F)


###############################################################################################################
####                              Plotting density of mutations with recombinations                        ####
###############################################################################################################

require(ggplot2)

window_recom_table = read.delim(output2, sep = "\t", header = T)

xbreaks = seq(0,chr_length,200000)
xlabels = seq(0,chr_length/1000,200000/1000)

size_dot = 0.2; size_axis_lines = 0.3; text_y_offset = 4; font = "Times"; dot_color = "dimgray";
axis_text_size = 10; axis_title_size = 15; ann_text_size = 5;
g1 <- ggplot(window_recom_table, aes(x=window_start, y=mean_density)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
  geom_line(size = size_dot, linetype = "solid") +
  scale_x_continuous(breaks = xbreaks, labels = xlabels) +
  # xlim(-1, chr_length) + 
  ylim(0, 20) + 
  xlab("Position (kilobases)") +
  ylab("Mean number of variants per kilobase") +
  ggtitle("Recombination density") +
  theme(text = element_text(family = font)) +
  theme(axis.text = element_text(size=axis_text_size, color="black"), axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))
plot_file = "recombination_density.pdf";
ggsave(plot_file, plot = g1, device = "pdf",width = 10, height = 5, dpi = 300, units = "in")




