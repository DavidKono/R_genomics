library(data.table)
library(dplyr)
library(doSNOW)
library(foreach)
library(parallel)

cleaned_flattened <- fread("ESC_complete_seqmonk.mapped.cleaned.flattened.bed", header = FALSE, sep = "\t")
cleaned_flattened <- na.omit(cleaned_flattened)
split_cleaned_flattened <- split(cleaned_flattened, cleaned_flattened$V1)


#Coordinates.Main_Chromosomes.All_Transcripts
split_coordinates <- split(Coordinates.Main_Chromosomes.All_Transcripts, Coordinates.Main_Chromosomes.All_Transcripts$chrom_name)


#initialise output to be n + 2 copy of cleaned_flattened
#chrA, startA, endA, chrB, startB, endB, AnnotationA, AnnotationB, TSS_Overlap_A, TSS_Overlap_B
output <- cleaned_flattened
output[, ':=' (TSS_Overlap_A = "", TSS_Overlap_B = "")]
setnames(output, old = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"), new = c("chrA", "startA", "endA", "chrB", "startB", "endB", "AnnotationA", "AnnotationB"))

contains_tss <- function(tss, start, end){
  return((start <= tss & tss <= end))
}

#intialise cluster
cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)

#parallelise using foreach on each chr
chr_results <- foreach(chr = names(split_cleaned_flattened), .combine = list, .multicombine = TRUE, .packages = c('data.table', 'dplyr')) %dopar% {
  #get subset of queries and reference with same chr
  reference <- as.data.table(split_cleaned_flattened[[chr]])
  query <- as.data.table(split_coordinates[[chr]]) 
  
  #define output for this chr
  chr_output <- reference[, .(chrA = V1, startA = V2, endA = V3, chrB = V4, startB = V5, endB = V6, AnnotationA = V7, AnnotationB = V8, TSS_Overlap_A = "", TSS_Overlap_B = "")]
  
  for (i in 1:nrow(reference)){

    #check if this reference has overlap with any query
    #get overlap indices
    A_overlaps_i <- contains_tss(query$TSS, reference$V2[i], reference$V3[i])
    B_overlaps_i <- contains_tss(query$TSS, reference$V5[i], reference$V6[i])
    
    #get overlaps
    
    A_overlaps <- paste(query$ensembl_gene_id[A_overlaps_i], 
                        query$gene_symbol[A_overlaps_i], 
                        query$gene_biotype[A_overlaps_i], 
                        query$ensembl_transcript_id[A_overlaps_i],
                        sep = ":")
    
    B_overlaps <- paste(query$ensembl_gene_id[B_overlaps_i], 
                        query$gene_symbol[B_overlaps_i], 
                        query$gene_biotype[B_overlaps_i], 
                        query$ensembl_transcript_id[B_overlaps_i],
                        sep = ":")
  
    #if multiple overlaps per reference, combine as string same field
    #else "No Match"
    #ensembl_gene_id:gene_symbol:gene_biotype:ensembl_transcript_id
    chr_output$TSS_Overlap_A[i] <- ifelse(length(A_overlaps) == 0, "No Match", paste(A_overlaps, collapse = ","))
    chr_output$TSS_Overlap_B[i] <- ifelse(length(B_overlaps) == 0, "No Match", paste(B_overlaps, collapse = ","))
    
  }
  #need to return sub-output when done
  return(chr_output)
}

stopCluster(cl)

output <- rbindlist(chr_results)
fwrite(output, "coordinates_overlaps.txt", sep = "\t")
fwrite(output, "coordinates_overlaps.bed", sep = "\t")

print("finished")

cat("length of output", nrow(output))
cat("num of output A that are non-matches",sum(output$TSS_Overlap_A == "No Match"))
cat("num of output B that are non-matches",sum(output$TSS_Overlap_B == "No Match"))







