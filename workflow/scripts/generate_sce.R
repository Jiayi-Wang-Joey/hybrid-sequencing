suppressPackageStartupMessages({
    library(Matrix)
    library(data.table)
    library(scuttle)
    library(SingleCellExperiment)
    library(Biostrings)
})

mtx <-  readMM(args$mtx)
mtx <- t(mtx)
bcd <- fread(args$bcd, header=FALSE)
fts <- fread(args$fts, header=FALSE)
if (TRUE) {
    i <- substr(wcs$sample, 1, 7)
    p <- paste0("/home/jiayiwang/loggedfs_david_penton/result/cellranger/",
                        i, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
    whitelist <- fread(p, header = FALSE)$V1
    whitelist <- as.character(
        reverseComplement(DNAStringSet(sub("-1$", "", whitelist))))
    
}
keep_bc <- bcd$V1 %in% whitelist
mtx <- mtx[, keep_bc]
bcd <- bcd[keep_bc, , drop = FALSE]
sce <- SingleCellExperiment(
    assays=list(counts=mtx),
    colData=DataFrame(sample=rep(wcs$sample, nrow(bcd)),
                      condition=rep(gsub("[0-9]", "", wcs$sample), nrow(bcd)),
                      patient=rep(gsub("[^0-9]", "", wcs$sample), nrow(bcd))))
rownames(sce) <- fts$V1
colnames(sce) <- bcd$V1

gene_mapping <- vroom::vroom(args$tx)
ids <- gene_mapping$gene_id[match(rownames(sce), gene_mapping$transcript_id)]

gsce <- aggregateAcrossFeatures(sce, ids = ids)
#gsce <- gsce[, colSums(counts(gsce)) > 0]

saveRDS(sce, args$tsce)
saveRDS(gsce, args$gsce)

