library(data.table)
library(dplyr)

#read in files
test_regions <- fread("Test_regions.bed", header = FALSE, sep = "\t")

methylation_calls <- fread("Methylation_Calls.Subsampling_100000_calls.E4.5-5.5.bed", header = TRUE, sep = "\t")
setnames(methylation_calls, old = "#chr", new = "chr")


#get overlaps, then match chromatin numbers
setkey(test_regions, V2, V3)
setkey(methylation_calls, start, end)

overlaps <- foverlaps(test_regions, methylation_calls, nomatch = 0)

chr_overlaps <- overlaps[chr == V1]

chr_overlaps[, Region := paste(chr_overlaps$V1, chr_overlaps$V2, chr_overlaps$V3, sep = ",")]


#initialise met_reads
region_names <- vector("character", nrow(test_regions))
for (i in 1:nrow(test_regions)) {
  Region <- test_regions[i, ]
  region_names[i] <- paste(Region$V1, Region$V2, Region$V3, sep = ",")
}
met_reads <- data.table(Region = region_names)
cell_names <- unique(methylation_calls$cell)
for (cell_name in cell_names) {
  met_reads[, (cell_name) := 0]
}

non_met_reads <- copy(met_reads)


#add met read count for each region at the specified cell site
for (cell_name in unique(chr_overlaps$cell)) {
  #match cell from overlaps to cell name, then add met read counts from overlaps, matched on the region row
  met_reads[chr_overlaps[cell == cell_name], (cell_name) := get(cell_name) + i.met_reads, on = .(Region)]
  non_met_reads[chr_overlaps[cell == cell_name], (cell_name) := get(cell_name) + i.non_met_reads, on = .(Region)]
}


fwrite(met_reads, "Cellular_Matrix_Met_Reads.bed")
fwrite(non_met_reads, "Cellular_Matrix_Non_Met_Reads.bed")


print("Finished")





