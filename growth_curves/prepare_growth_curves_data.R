
# This R script is used to prepare data to plot bacterial growth curves using R

rawdata <- read.csv("curves_all_14-10-22.csv")
dim(rawdata)
# [1]  48 598

### R functions to extract sample id, repeat id, and daptomycin concentration from raw data table
extract_sample_id = function(x)
{
  # "X10770_3.44_000.2"
  # "MPROS2026_000.3"
  sample_id = "";
  tmp = unlist(strsplit(x, "\\_"));
  sample_id = paste(tmp[-length(tmp)], collapse = "_");
  return(sample_id)
}

extract_repeat_id = function(x)
{
  repeat_id = 0; 
  if(grepl("\\.[[:digit:]]$", x))
  {
    tmp = unlist(strsplit(x, "\\."));
    repeat_id = as.numeric(tmp[length(tmp)]) + 1; 
  } else { repeat_id = 1; }
  return(repeat_id)
}

extract_environment_id = function(x)
{
  environment = "";
  tmp = unlist(strsplit(x, "\\_"));
  environment = unlist(strsplit(tmp[length(tmp)], "\\."))[1];
  return(environment)
}

metadata = mat.or.vec(ncol(rawdata)-1, 4)
colnames(metadata) = c("curve_id","sample_id","repeat_id","environment")
metadata[,1] = as.vector(colnames(rawdata)[-1])
metadata[,2] = as.vector(sapply(colnames(rawdata)[-1], extract_sample_id))
metadata[,3] = as.vector(sapply(colnames(rawdata)[-1], extract_repeat_id))
metadata[,4] = as.vector(sapply(colnames(rawdata)[-1], extract_environment_id))
# fixing sample ids
metadata[,2] = gsub("^X", "", metadata[,2])
metadata[,2] = gsub("\\.", "\\#", metadata[,2])

write.table(metadata, file = "curves_all_14-10-22.curves_metadata.csv", sep = "\t", col.names = T, row.names = F, quote = F)

unique_sample_ids = unique(as.vector(metadata[,2]))
unique_sample_ids[which(unique_sample_ids=="")] = "negative_control"

write.table(unique_sample_ids, file = "curves_all_14-10-22.strain_metadata.csv", sep = "\n", col.names = F, row.names = F, quote = F)
# NOTE: the file curves_all_14-10-22.strain_metadata.csv was manually filled in 

