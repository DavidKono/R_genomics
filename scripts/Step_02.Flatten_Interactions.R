cleaned <- read.table("ESC_complete_seqmonk.mapped.cleaned.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

flattened_output <- file("ESC_complete_seqmonk.mapped.cleaned.flattened.bed", "w")

for (i in seq(1, nrow(cleaned), by = 2)) {
  first_row <- cleaned[i,]
  second_row <- cleaned[i+1,]
  
  flattened <- c(first_row[1:3], second_row[1:3], first_row[4], second_row[4])
  flattened_string <- paste(flattened, collapse = "\t")    #tab separation
  writeLines(flattened_string, flattened_output)
}

close(flattened_output)
cleaned_flattened <- read.table("ESC_complete_seqmonk.mapped.cleaned.flattened.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

print("Finished flattening")

