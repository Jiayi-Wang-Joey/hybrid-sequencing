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

sce <- SingleCellExperiment(
    assays=list(counts=mtx),
    colData=DataFrame(sample=rep(wcs$sample, nrow(bcd)),
                      condition=rep(gsub("[0-9]", "", wcs$sample), nrow(bcd)),
                      patient=rep(gsub("[^0-9]", "", wcs$sample), nrow(bcd))))
rownames(sce) <- fts$V1
colnames(sce) <- bcd$V1
counts(sce) <- as(counts(sce), "dgCMatrix")
gene_mapping <- vroom::vroom(args$tx)
ids <- gene_mapping$gene_id[match(rownames(sce), gene_mapping$transcript_id)]

#gsce <- aggregateAcrossFeatures(sce, ids = ids)
#gsce <- gsce[, colSums(counts(gsce)) > 0]

saveRDS(sce, args$tsce)
#saveRDS(gsce, args$gsce)