library(data.table)
library(dplyr)

test_regions <- fread("Test_regions.bed", header = FALSE, sep = "\t")
#test_regions <- fread("Test_regions.Subsample_100.bed", header = FALSE, sep = "\t")
setnames(test_regions, c("V1", "V2", "V3"), c("chr", "start", "end"))

methylation_calls <- fread("Methylation_Calls.Subsampling_100000_calls.E4.5-5.5.bed", header = TRUE, sep = "\t")
split_methylation_calls <- split(methylation_calls, methylation_calls$cell)


#initialise met_reads
region_names <- vector("character", nrow(test_regions))
for (i in 1:nrow(test_regions)) {
  region <- test_regions[i, ]
  region_names[i] <- paste(region$chr, region$start, region$end, sep = ",")
}
met_reads <- data.table(region = region_names)
cell_names <- names(split_methylation_calls)
for (cell_name in cell_names) {
  met_reads[, (cell_name) := 0]
}


names(split_methylation_calls)

setkey(test_regions, start, end)
for (cell_name in cell_names){
  cell <- split_methylation_calls[[cell_name]]
  
  
  #sewt key region$V2, region$V3, cell$start, cell$end)
  setkey(cell, start, end)
  overlaps <- foverlaps(test_regions, cell, nomatch = 0)
  
  overlaps_sum <- overlaps[, .(total_value = sum(met_reads)), by = .(region)]
  
  #print(met_reads[[cell_name]][i])
    

  print(cell_name)
  
  print(sum(cell$cell_name))
  #test first cell only
  Sys.sleep(5)
  break
}
