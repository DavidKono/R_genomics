#ensure first get formatted bed file from process_data.R

#install.packages("data.table")
library(data.table)
#install.packages("dplyr")
library(dplyr)

cleaned_flattened <- fread("ESC_complete_seqmonk.mapped.cleaned.flattened.bed", header = FALSE, sep = "\t")
cleaned_flattened <- na.omit(cleaned_flattened)
split_cleaned_flattened <- split(cleaned_flattened, cleaned_flattened$V1)

#query_file <- "scATAC_formatted.bed"
query_file <- "scChipSeq_formatted.bed"
query_formatted <- fread(query_file, header = FALSE, sep = "\t")
query_formatted <- na.omit(query_formatted)
split_query_formatted <- split(query_formatted, query_formatted$V1)


#create output file based on input bed
file_prefix <- sub("\\.bed$", "", query_file)
output_filename <- paste0(file_prefix, "_output.bed")

#mutate original for output
output <- query_formatted %>%
  mutate(RegionC = paste(V1, V2, V3, sep = ","),
      RegionD = paste(V4, V5, V6, sep = ","),
      AnnotationsQ = V7,
      RegionA = "No Full Match",
      RegionB = "No Full Match",
      AnnotationsR = "No Full Match") %>%
  select(RegionC, RegionD, AnnotationsQ, RegionA, RegionB, AnnotationsR)

query_overlaps <- data.table()
total_overlaps <- data.table()

for (chr in names(split_cleaned_flattened)) {
  print(chr)
  reference <- as.data.table(split_cleaned_flattened[[chr]])
  query <- as.data.table(split_query_formatted[[chr]])
  
  #region 1 on region 1 ie A and C
  setkey(reference, V2, V3)
  setkey(query, V2, V3)
  region1_on1 <- foverlaps(reference, query, nomatch = 0)
  
  #region 1 on region 2 ie A and D
  setkey(reference, V2, V3)
  setkey(query, V5, V6)
  region1_on2 <- foverlaps(reference, query, nomatch = 0)
    
  #region 2 on region 2 ie B and D
  setkey(reference, V5, V6)
  setkey(query, V5, V6)
  region2_on2 <- foverlaps(reference, query, nomatch = 0)
  
  #region 2 on region 1 ie B and D
  setkey(reference, V5, V6)
  setkey(query, V2, V3)
  region2_on1 <- foverlaps(reference, query, nomatch = 0)
  
  Aoverlaps <- rbind(region1_on1, region1_on2)
  Boverlaps <- rbind(region2_on1, region2_on2)
  
  all_overlaps <- rbind(Aoverlaps, Boverlaps)

  #get where both regions have overlaps
  overlaps <- merge(Aoverlaps, Boverlaps)

  #first case of query overlapping with reference
  current_query_overlaps <- overlaps %>%
    group_by(V1, V2, V3, V4, V5, V6, V7) %>%  #take query
    slice(1) %>%     #take first occurance
    ungroup() 
  
  current_overlaps <- all_overlaps %>%
    group_by(V1, V2, V3, V4, V5, V6, V7) %>%  #take query
    slice(1) %>%     #take first occurance
    ungroup() 

  total_overlaps <- bind_rows(total_overlaps, current_overlaps)
  
  query_overlaps <- bind_rows(query_overlaps, current_query_overlaps)
  
  print(nrow(query_overlaps))
}

overlapped_formatted <- query_overlaps %>%
  mutate(RegionC = paste(V1, V2, V3, sep = ","),
    RegionD = paste(V4, V5, V6, sep = ","),
    AnnotationsQ = V7,
    RegionA = paste(i.V1, i.V2, i.V3, sep = ","),
    RegionB = paste(i.V4, i.V5, i.V6, sep = ","),
    AnnotationsR = paste(i.V7, V8, sep = ",")) %>%
  select(RegionC, RegionD, AnnotationsQ, RegionA, RegionB, AnnotationsR)

output <- rows_update(output, overlapped_formatted, by = c("RegionC", "RegionD", "AnnotationsQ"))
fwrite(output, output_filename, append = TRUE, sep = "\t", col.names = !file.exists(output_filename))

#stats
cat("Query file being checked:", query_file)

cat("Interactions (pairs) in Reference:", nrow(cleaned_flattened), "\n")
cat("Interactions (pairs) in Query:", nrow(query_formatted), "\n")

cat("Query with match:", nrow(query_overlaps), "\n")
cat("Query without match:", nrow(query_formatted) - nrow(query_overlaps), "\n")


print("Finished getting overlaps")


