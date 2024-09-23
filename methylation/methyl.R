library(data.table)
library(dplyr)

test_regions <- fread("Test_regions.bed", header = FALSE, sep = "\t")
#test_regions <- fread("Test_regions.Subsample_100.bed", header = FALSE, sep = "\t")
test_regions[, position := .I]

methylation_calls <- fread("Methylation_Calls.Subsampling_100000_calls.E4.5-5.5.bed", header = TRUE, sep = "\t")
setnames(methylation_calls, old = "#chr", new = "chr")
methylation_calls[, position := .I]
split_methylation_calls <- split(methylation_calls, methylation_calls$cell)


#initialise met_reads
region_names <- vector("character", nrow(test_regions))
for (i in 1:nrow(test_regions)) {
  region <- test_regions[i, ]
  region_names[i] <- paste(region$V1, region$V2, region$V3, sep = ",")
}
met_reads <- data.table(region = region_names)
cell_names <- names(split_methylation_calls)
for (cell_name in cell_names) {
  met_reads[, (cell_name) := 0]
}

overlap <- function(start1, end1, start2, end2) {
  return((start1 == start2 | end1 == start2 | end1 == end2 | start1 == end2))
}


#first split among cells
for (cell_name in cell_names){
    cell <- split_methylation_calls[[cell_name]]
    
    #next split among chr
    cell_chr_split <- split(cell, cell$chr)
    test_regions_split <- split(test_regions, test_regions$V1)
    
    for (chr_name in names(cell_chr_split)){
      cell_chr <- cell_chr_split[[chr_name]]
      test_regions_chr <- test_regions_split[[chr_name]]
      
      for (i in 1:nrow(test_regions_chr)) {
        region <- test_regions_chr[i]
        
        #see which have overlaps
        overlap_indices <- overlap(region$V2, region$V3, cell_chr$start, cell_chr$end)
        #get actual position of any overlaps
        actual_positions <- cell_chr$position
        true_indices <- actual_positions[which(overlap_indices)]
        
        overlaps <- sum(cell$met_reads[overlap_indices])
        
        #print(overlaps)
        #print(cell_name)
        actual_indice <- region$position
        met_reads[[cell_name]][actual_indice] <- overlaps
        
        #print(i)
        #print(met_reads[[cell_name]][i]) 
      }
      
    }
    
  print(sum(met_reads$cell_name))
}


print("finished")







