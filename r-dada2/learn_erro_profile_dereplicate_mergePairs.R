
learn_error_dereplicate_merge_pairs <- function(path) {
  filtFs <- sort(list.files(path, pattern="_F_filt.fastq.gz", full.names = TRUE))
  filtRs <- sort(list.files(path, pattern="_R_filt.fastq.gz", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
  # Learn the error rates
  set.seed(100)
  errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 15)
  errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 15)
  saveRDS(errF, file = paste(path, "errF_", unlist(strsplit(path, "/"))[3], ".rds", sep = ""))
  saveRDS(errR, file = paste(path, "errR_", unlist(strsplit(path, "/"))[3], ".rds", sep = ""))
  print(paste("Done with learnErrors for:", path))
  #plotErrors(errF, nominalQ=TRUE)
  # Dereplecation
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  print(paste("Done with dereplication for:", path))
  # Name the derep-class objects by the sample names
  # names(derepFs) <- sample.names
  # names(derepRs) <- sample.names
  # Sample Inference
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  saveRDS(dadaFs, file = paste(path, "dadaFs_", unlist(strsplit(path, "/"))[3], ".rds", sep = ""))
  saveRDS(dadaRs, file = paste(path, "dadaRs_", unlist(strsplit(path, "/"))[3], ".rds", sep = ""))
  #dadaFs
  print(paste("Done with sample inference for:", path))
  # Merge paired reads
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  print(paste("Done with pair-end merger for:", path))
  # Inspect the merger data.frame from the first sample
  #head(mergers[[1]])
  # Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, file = paste(path, "seqtab_", unlist(strsplit(path, "/"))[3], ".rds", sep = ""))
  dim(seqtab)
  table(nchar(getSequences(seqtab)))
  print(paste("Done with:", path))
}


