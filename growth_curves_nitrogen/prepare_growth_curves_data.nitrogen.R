
# This R script is used to prepare data to plot bacterial growth curves using R for nitrogen biolog plates
# Output *.curves_metadata.csv files must contain: "curve_id","sample_id","repeat_id" and "environment"
# NOTE: sample 103356 lacks repeat id in Export_All_Plates_Single_File.csv. Manually added in Export_All_Plates_Single_File.edited.csv

rawdata <- read.csv("Export_All_Plates_Single_File.edited.csv", sep = ",", header = T)
dim(rawdata)
# [1] 5819  105

# removing duplicated header
tmp = which(rawdata$Plate.File == "Plate File")
rawdata = rawdata[-tmp,]
dim(rawdata)
# [1] 5790  105


# Exploring and understanding the data
sort(table(rawdata$Field.2))
# 103356 A 103356 B 103356 C 144398 A 144398 B 144398 C 146083 A 146083 B 146083 C  58870 A  58870 B 
#      193      193      193      193      193      193      193      193      193      193      193 
# 58870 C  74975 A  74975 B  74975 C  80131 A  80131 B  80131 C    JE2 A    JE2 B    JE2 C NE1113 A 
#     193      193      193      193      193      193      193      193      193      193      193 
# NE1113 B NE1113 C  NE286 A  NE286 B  NE286 C  NE857 A  NE857 B  NE857 C 
#     193      193      193      193      193      193      193      193 


# NOTE: Field 2 contains strain id and replicate (A, B, and C)
# NOTE: columns A01 to H12 include the different nitrogen substrate

environment_ids = colnames(rawdata)[10:ncol(rawdata)]

extract_sample_id = function(x)
{
  # 144398 A
  sample_id = "";
  tmp = unlist(strsplit(x, " "));
  sample_id = tmp[1]
  return(sample_id)
}

extract_repeat_id = function(x)
{
  # 144398 A
  repeat_id = "";
  tmp = unlist(strsplit(x, " "));
  repeat_id = tmp[2]
  return(repeat_id)
}

rawdata$sample_id = sapply(rawdata$Field.2, extract_sample_id)
rawdata$repeat_id = sapply(rawdata$Field.2, extract_repeat_id)

table(rawdata$sample_id)
# 103356 144398 146083  58870  74975  80131    JE2 NE1113  NE286  NE857 
#    579    579    579    579    579    579    579    579    579    579 

table(rawdata$repeat_id)
#    A    B    C 
# 1930 1930 1930 


## Creating curve ids, curve metadata and re-formatting raw curve data
sample_ids = sort(unique(rawdata$sample_id))
repeat_ids = sort(unique(rawdata$repeat_id))
hours = sort(unique(as.numeric(rawdata$Hour)))
num_rows = length(sample_ids)*length(repeat_ids)*length(environment_ids)
# [1] 2880

metadata = mat.or.vec(1, 4)
growthdata = hours; # time points as rows, curves as columns
colnames(metadata) = c("curve_id","sample_id","repeat_id","environment")
curve_ids = vector()

for(s in 1:length(sample_ids)){
  for(r in 1:length(repeat_ids)){
    for(e in 1:length(environment_ids)){
      sample_id = sample_ids[s]
      repeat_id = repeat_ids[r]
      environment_id = environment_ids[e]
      # creating curve metadata entry
      curve_id = paste(sample_id, repeat_id, environment_id, sep = "_")
      curve_ids = c(curve_ids, curve_id)
      metadata = rbind(metadata, c(curve_id, sample_id, repeat_id, environment_id))
      # extracting measured growth data for curve_id
      r_idx = which(rawdata$sample_id == sample_id & rawdata$repeat_id == repeat_id)
      rawdata_curve = as.vector(rawdata[r_idx,environment_id])
      growthdata = cbind(growthdata, rawdata_curve)
    }
  }
}
metadata = metadata[-1,]
dim(metadata)
# [1] 2880    4
length(curve_ids)
# [1] 2880
dim(growthdata)
# [1]  193 2881
colnames(growthdata) = c("time", curve_ids)


## Adding nitrogen source information
# input file created manually from 00A-042-Rev-C-Phenotype-MicroArrays-1-10-Plate-Maps.pdf
pm3 = read.delim("00A-042-Rev-C-Phenotype-MicroArrays-1-10-Plate-Maps.PM3.csv", sep = " ", header = F)
colnames(pm3) = c("environment", "nitrogen_source_pm3")
# making nitrogen source match
code = metadata[,"environment"]
code = gsub("A0", "A", code)
code = gsub("B0", "B", code)
code = gsub("C0", "C", code)
code = gsub("D0", "D", code)
code = gsub("E0", "E", code)
code = gsub("F0", "F", code)
code = gsub("G0", "G", code)
code = gsub("H0", "H", code)
metadata = cbind(metadata, code)
metadata2 = merge(metadata, pm3, by.x = "code", by.y = "environment", all.x = TRUE)
dim(metadata2)
# [1] 2880    5

metadata2 = metadata2[order(metadata2$curve_id), ]

write.table(metadata2, file = "Export_All_Plates_Single_File.curves_metadata.csv", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(growthdata, file = "Export_All_Plates_Single_File.curves_growthdata.csv", sep = "\t", col.names = T, row.names = F, quote = F)



# to do

write.table(sample_ids, file = "Export_All_Plates_Single_File.strain_metadata.csv", sep = "\n", col.names = F, row.names = F, quote = F)
# NOTE: the file curves_all_14-10-22.strain_metadata.csv was manually filled in 

