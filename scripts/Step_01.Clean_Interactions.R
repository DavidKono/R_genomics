mapped <- read.table("ESC_complete_seqmonk.mapped.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
unmapped <- read.table("ESC_complete_seqmonk.unmapped.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

cleaned_output <- file("ESC_complete_seqmonk.mapped.cleaned.bed", "w")

pair_array = array(dim = nrow(unmapped))
for (i in 1:nrow(unmapped)) {
  v4 <- unmapped[i,4]
  #print(v4)
  fourth <- strsplit(v4, ":")
  pair <- sapply(fourth, function(x) x[length(x)])
  #print(pair)
  pair_array[i] <- pair
}
print(pair_array)


for (i in 1:nrow(mapped)) {
  v4 <- mapped[i,4]
  fourth <- strsplit(v4, ":")
  pair <- sapply(fourth, function(x) x[length(x)])
  
  if (!(pair %in% pair_array)) { 
    row_string <- paste(mapped[i, ], collapse = "\t")   #need as string
    writeLines(row_string, cleaned_output)
  }
}

cleaned <- read.table("ESC_complete_seqmonk.mapped.cleaned.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

print("Finished cleaning")

