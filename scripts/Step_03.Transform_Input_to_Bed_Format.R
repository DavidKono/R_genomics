library(data.table)
table <- read.table("scATAC_filtered_connections_cicero_0.25.txt",header = TRUE,stringsAsFactors=FALSE, quote="")
#table <- read.table("scChipSeq_filtered_connections_cicero_0.25.txt",header = TRUE,stringsAsFactors=FALSE, quote="")

output_file <- "scATAC_formatted.bed"
#output_file <- "scChipSeq_formatted.bed"
formatted_lines <- data.table()

for (i in 1:nrow(table)) {
  row <- table[i, ]
  
  chr1_str <- gsub('"', '', row[1])
  chr2_str <- gsub('"', '', row[2])
  value <- as.numeric(row[3])
  
  #split further to get chr number, start, end
  chr1_parts <- strsplit(chr1_str, "_")[[1]]
  chr2_parts <- strsplit(chr2_str, "_")[[1]]
  
  formatted_line <- data.table(
    chr1_parts[1], as.numeric(chr1_parts[2]), as.numeric(chr1_parts[3]),
    chr2_parts[1], as.numeric(chr2_parts[2]), as.numeric(chr2_parts[3]),
    value
  )
  
  formatted_lines <- rbind(formatted_lines, formatted_line)
}

fwrite(formatted_lines, output_file, sep = "\t", col.names = FALSE)

print("Finished reformatting")