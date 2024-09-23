library(data.table)
library(dplyr)

methylation_calls <- fread("Methylation_Calls.Subsampling_100000_calls.E4.5-5.5.bed", header = TRUE, sep = "\t")
setnames(methylation_calls, old = "#chr", new = "chr")

REPTILE <- fread("REPTILE_mESC_Training_Regions.bed", header = FALSE, sep = "\t")
REPTILE <- data.table(REPTILE)

#ensure 2k bps long
for (i in 1:nrow(REPTILE)){
  center <- round((REPTILE$V2[i] + REPTILE$V3[i]) / 2 )
  REPTILE$V2[i] <- center - 1000  
  REPTILE$V3[i] <- center + 1000  
}

# test data
# set.seed(1)
# 
# num_rows <- 100
# 
# chromosomes <- sample(1:20, num_rows, replace = TRUE)
# chromosomes <- paste0("chr", chromosomes)
# 
# start_sites <- sample(1:(1e7 - 2000), num_rows)
# create regions in the format chr1_100000_102x000
# regions <- paste0(chromosomes, "_", start_sites, "_", start_sites + 2000)
# 
# Candidate_Enhancer_Regions <- data.table(
#   Region = regions
# )

#create name ids
Candidate_Enhancer_Regions <- copy(REPTILE)
Candidate_Enhancer_Regions[, Region := paste0(V1, "_", V2, "_", V3, "_", V4)]
Candidate_Enhancer_Regions[, V1 := NULL]
Candidate_Enhancer_Regions[, V2 := NULL]
Candidate_Enhancer_Regions[, V3 := NULL]
Candidate_Enhancer_Regions[, V4 := NULL]

Candidate_Met_Counts <- copy(Candidate_Enhancer_Regions)
Candidate_Non_Met_Counts <- copy(Candidate_Enhancer_Regions)

#add window headers
window_names <- vector("character", 20)
met_count_names <- vector("character", 20)
non_met_count_names <- vector("character", 20)

for (i in 1:20) {
  num <- i
  window_names[i] <- paste("Window", num, "Start", sep = "_")
  met_count_names[i] <- paste("Window", num, "Met_Count", sep = "_")
  non_met_count_names[i] <- paste("Window", num, "Non_Met_Count", sep = "_")
}


for (i in 1:20) {
  Candidate_Enhancer_Regions[, (window_names[i]) := 0]
  Candidate_Met_Counts[, (met_count_names[i]) := 0]
  Candidate_Non_Met_Counts[, (non_met_count_names[i]) := 0]
}
Candidate_Enhancer_Regions[, ("Window_End") := 0]


#20 loop of taking last collumn, adding 100
for (i in 1:nrow(Candidate_Enhancer_Regions)){
  Candidate_Enhancer_Regions_split <- strsplit(Candidate_Enhancer_Regions[[1]][i], "_")
  window_start = Candidate_Enhancer_Regions_split[[1]][2]
  window_start <- as.numeric(window_start)
  
  for (j in 1:21){
    Candidate_Enhancer_Regions[[j + 1]][i] <- window_start
    window_start <- window_start + 100
  }
}

#only chr19 so dont subdivide for now
chrs <- unique(REPTILE$V1)
for (chr_num in chrs){
  for (i in 1:20) {
    methylation_calls_chr <- methylation_calls[chr == chr_num]
    
    setkey(methylation_calls_chr, start, end)
    collumns <- c(names(Candidate_Enhancer_Regions)[i + 1], names(Candidate_Enhancer_Regions)[i + 2])
    print(collumns)
    setkeyv(Candidate_Enhancer_Regions, collumns)
    overlaps <- foverlaps(methylation_calls_chr, Candidate_Enhancer_Regions, nomatch = 0)
    
    #add met_counts where same id is in data table
    Candidate_Met_Counts[overlaps, on = .(Region), (met_count_names[i]) := get(met_count_names[i]) + i.met_reads]
    Candidate_Non_Met_Counts[overlaps, on = .(Region), (non_met_count_names[i]) := get(non_met_count_names[i]) + i.non_met_reads]
  }
}



#format Candidate_Enhancer_Regions again for output ie add chr num
#complete resuse of code not ideal, but better than adding chr num before processing
for (i in 1:nrow(Candidate_Enhancer_Regions)){
  Candidate_Enhancer_Regions_split <- strsplit(Candidate_Enhancer_Regions[[1]][i], "_")
  
  chr_num <- Candidate_Enhancer_Regions_split[[1]][1]
  
  window_start <- Candidate_Enhancer_Regions_split[[1]][2]
  window_start <- as.numeric(window_start)
  
  for (j in 1:21){
    Candidate_Enhancer_Regions[[j + 1]][i] <- paste0(chr_num, ":", window_start)
    window_start <- window_start + 100
  }
}

fwrite(Candidate_Enhancer_Regions, "Candidate_Enhancer_Regions.bed", sep = "\t", col.names = TRUE)

fwrite(Candidate_Met_Counts, "Met_Reads_Output.bed", sep = "\t", col.names = TRUE)
fwrite(Candidate_Non_Met_Counts, "Non_Met_Reads_Output.bed", sep = "\t", col.names = TRUE)
fwrite(Candidate_Met_Counts, "Met_Reads_Output.txt", sep = "\t", col.names = TRUE)
fwrite(Candidate_Non_Met_Counts, "Non_Met_Reads_Output.txt", sep = "\t", col.names = TRUE)

print("Finished")


