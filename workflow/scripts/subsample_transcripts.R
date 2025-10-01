suppressPackageStartupMessages(library(rtracklayer))
transcriptome <- import(args$gtf)
transcript_ids <- unique(transcriptome$transcript_id)
sampled_transcripts <- sample(x = transcript_ids, 
                              replace = FALSE, 
                              size = args[[1]])
subsampled_transcriptome <- transcriptome[transcriptome$transcript_id %in% sampled_transcripts]
export(
    subsampled_transcriptome,
    args$res
)
