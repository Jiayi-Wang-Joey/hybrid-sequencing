suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicAlignments)
    library(dplyr)
    library(biomaRt)
})

homopolymerFinder <- function(seq)
    # Finds homopolymers in a given DNAStringSet object.
{
    if (!is(seq, "DNAStringSet")) {
        stop("'seq' should be a DNAStringSet object")
    }
    
    all.homo <- .Call(cxx_find_homopolymers, seq)
    homo.range <- IRanges(all.homo[[2]], all.homo[[2]] + all.homo[[3]] - 1L)
    mcols(homo.range)$base <- all.homo[[4]]
    
    groupings <- factor(all.homo[[1]]+1L, levels=seq_along(seq))
    output <- split(homo.range, groupings, drop=FALSE)
    names(output) <- names(seq)
    return(output)
}

.convert_to_single_letter <- function(protein_change, save_p=FALSE) {
    aa_mapping <- c(
        Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", 
        Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I", 
        Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P", 
        Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V"
    )
    
    if (!save_p) protein_change <- gsub("p\\.", "", protein_change)
    for (aa in names(aa_mapping)) {
        protein_change <- gsub(aa, aa_mapping[[aa]], protein_change)
    }
    protein_change
}

# .extractAnnotation <- \(vcfs, id = "sample") {
#     annos <- lapply(vcfs, \(vcf) {
#         ann <- info(vcf)$ANN
#         df <- lapply(ann, \(a) {
#             ls <- lapply(a, \(x) {
#                 y <- strsplit(x, "\\|")[[1]]
#                 if ("missense_variant" %in% y) {
#                     res <- y[1:11]
#                 } else {
#                     res <- NULL
#                 }
#             })
#             ls <- ls[sapply(ls, \(x) !is.null(x))]
#             do.call(rbind, ls)
#         })
#         df <- data.frame(do.call(rbind,df))
#         names(df) <- c("Reference", "Variant_Type", "Impact",
#                        "Gene", "Gene_ID", "Feature_Type",
#                        "Transcript_ID", "Consequence",
#                        "Exon_Rank", "cDNA_Change",
#                        "Protein_Change")
#         data.table(df)
#     })
#     annos <- rbindlist(annos,idcol=id)
#     annos[,mutation:=paste0(Gene,":", .convert_to_single_letter(Protein_Change))]
#     annos[,mutation.p:=paste0(Gene,":",
#                               .convert_to_single_letter(Protein_Change, TRUE))]
#     annos[,nMutation:=uniqueN(mutation),by=sample]
#     annos[,condition:=substr(sample, 1, nchar(sample) - 1)]
#     annos[,patient:=paste0("patient",
#                            substr(sample, nchar(sample), nchar(sample)))]
#     annos
# }

.extractAnnotation <- \(vcf) {
    ann <- info(vcf)$ANN
    df <- lapply(ann, \(a) {
        ls <- lapply(a, \(x) {
            y <- strsplit(x, "\\|")[[1]]
            # if (any(grepl("missense_variant", y))) {
            #     res <- y[1:11]
            # } else {
            #     res <- NULL
            # }
        })
        ls <- ls[sapply(ls, \(x) !is.null(x))]
        if (length(ls)==0) print(i) 
        do.call(rbind, ls)
    })
    df <- lapply(df, \(.) .[1,])
    df <- do.call(rbind, df)
    colnames(df) <- c("Reference", "Variant_Type", "Impact",
                    "Gene", "Gene_ID", "Feature_Type",
                     "Transcript_ID", "Consequence",
                     "Exon_Rank", "cDNA_Change",
                      "Protein_Change")
    dt <- data.table(df)
    dt[,mutation:=paste0(Gene,":", .convert_to_single_letter(Protein_Change))]
    dt
}

.extractCoordinates <- function(input) {
    chr <- sub("^(.*):.*$", "\\1", input)                  
    gDNA <- sub(".*g\\.(\\d+)_?(\\d+)?/.*/.*", "\\1_\\2", input) 
    cDNA <- sub(".*c\\.(\\d+)_?(\\d+)?/.*", "\\1_\\2", input) 
    
    gDNA_coords <- unlist(strsplit(gDNA, "_"))
    cDNA_coords <- unlist(strsplit(cDNA, "_"))
    
    gDNA_start <- as.numeric(gDNA_coords[1])
    gDNA_end <- as.numeric(gDNA_coords[2])
    cDNA_start <- as.numeric(cDNA_coords[1])
    cDNA_end <- as.numeric(cDNA_coords[2])
    
    data.table(
        chr = chr,
        gDNA_start = gDNA_start,
        gDNA_end = gDNA_end,
        cDNA_start = cDNA_start,
        cDNA_end = cDNA_end
    )
}


.is_missense <- \(vcf) {
    ann_field <- info(vcf)$ANN
    if (nrow(vcf) != 0) {
        is_missense_variant <- \(ann) {
            any(grepl("missense_variant", ann, ignore.case = TRUE))}
        vcf <- vcf[sapply(ann_field, is_missense_variant)]
    }
    
    if(nrow(vcf)==0) {
        NULL
    } else {
        vcf
    }
}

findMutation <- \(bamFile, bamIndex, 
                  fastaFile, 
                  mutChr=NULL, 
                  mutStart=NULL, 
                  mutEnd=NULL,
                  gene, mutation,
                  translate = FALSE,
                  annotate = FALSE,
                  cellBarcode=NULL,
                  cellType) {
    # if (any(is.null(c(mutChr, mutStart, mutEnd)))) {
    #     aa <- paste0(gene, ":p.", mutation)
    #     aa <- substr(aa, 1, nchar(aa)-1)
    #     pos <- system(paste0("transvar panno --ccds -i", aa), intern=TRUE)[2][1]
    #     pos <- strsplit(pos, "\t")[[1]][5]
    #     pos <- .extractCoordinates(pos)
    #     mutChr <- pos$chr
    #     mutStart <- pos$gDNA_start
    #     mutEnd <- pos$gDNA_end
    # } 
    if (any(is.na(c(mutChr, mutStart, mutEnd)))) {
        return(NULL)
    }
    region <- GRanges(mutChr, IRanges(start=mutStart, end=mutEnd))
    fa <- FaFile(fastaFile)
    refSeq <- getSeq(fa, region)
    param <- ScanBamParam(which = region, 
                          what = ("seq"),
                          tag = "CB")
    ga <- readGAlignments(bamFile, bamIndex, param=param)
    mcols(ga)$mutSeq <- stackStringsFromGAlignments(ga,region)
    mcols(ga)$refSeq <- rep(as.character(refSeq[[1]]), length(ga))
    mcols(ga)$seq <- NULL
    ga <- as.data.table(ga)
    ga$cigar <- NULL
    ga <- ga[mutSeq != "..."]
    if (translate) {
        mutSeq <- sapply(ga$mutSeq, 
                         \(.) translate(DNAString(.)))
        ga$mutAA <- sapply(mutSeq, as.character)
        ga$refAA <- as.character(translate(refSeq[[1]]))
        ga[,ifMutate:=ifelse(mutAA==refAA, 0, 1)]
    } else {
        ga[,ifMutate:=ifelse(mutSeq==refSeq, 0, 1)]
    }
    ga[,gene:=gene]
    ga[,mutation:=mutation]
    ga[,VAF:=sum(ifMutate)/length(ifMutate)]
    ga[,coverage:=length(ifMutate)]
    if (annotate) {
        idx <- match(ga$CB, cellBarcode)
        ga$cellType <- as.character(cellType[idx])
        ga$cellType[is.na(ga$cellType)] <- "disgarded"
    }
    ga
}


findMutationSample <- \(bamFile, bamIndex,
                        mutStart, mutEnd, mutChr, 
                        refSeq, refMut) {
    region <- GRanges(mutChr, IRanges(start=mutStart, end=mutEnd))
    param <- ScanBamParam(which = region, 
                          what = ("seq"))
    ga <- readGAlignments(bamFile, bamIndex, param=param)
    mcols(ga)$mutSeq <- stackStringsFromGAlignments(ga,region)
    mcols(ga)$seq <- NULL
    ga <- as.data.table(ga)
    ga$cigar <- NULL
    ga <- ga[mutSeq != "."]
    ga[,read_depth:=.N]
    ga[,refSeq:=refSeq]
    ga[,nRef:=sum(mutSeq==refSeq)]
    ga[,nAlt:=sum(mutSeq!=refSeq)]
    ga[,VAF:=nAlt/read_depth]
    unique(ga[,.(nRef, nAlt, read_depth)])
}

filterByCell <- \(bamFile, bamIndex, 
                  start, end, chr, 
                  keep = c("Adrenocortical", "Subcapsular"),
                  cellBarcode=NULL,
                  cellType=NULL) {
    region <- GRanges(chr, IRanges(start=start, end=end))
    param <- ScanBamParam(which = region, 
                          what = c("seq"),
                          tag = "CB",
                          flag=scanBamFlag(isSecondaryAlignment = FALSE))
    ga <- readGAlignments(bamFile, bamIndex, param=param)
    ga <- as.data.table(ga)
    idx <- match(ga$CB, cellBarcode)
    cell <- as.character(cellType[idx])
    if (any(keep %in% cell)) {
        TRUE
    } else {
        FALSE
    }
}

.ifInsert <- function(cigar, start, end, target_pos) {
    # Split CIGAR string into operations
    ops <- unlist(strsplit(cigar, "(?<=[MIDNSHP])", perl = TRUE))
    
    current_pos <- start
    
    for (i in seq_along(ops)) {
        op <- ops[i]
        op_length <- as.numeric(substr(op, 1, nchar(op) - 1))
        op_type <- substr(op, nchar(op), nchar(op))
        
        if (op_type %in% c("M", "D", "N")) {
            if (current_pos <= target_pos && target_pos < current_pos + op_length) {
                # Check the next operation if it exists
                if (i + 1 <= length(ops) && grepl("I", ops[i + 1])) {
                    return(TRUE)
                }
            }
            current_pos <- current_pos + op_length
        } else if (op_type %in% c("I", "S", "H")) {
            next
        }
    }
    return(FALSE)
}

.ifDelete <- function(cigar, read_start, read_end, target_start, target_end) {
    # Split CIGAR string into operations
    ops <- unlist(strsplit(cigar, "(?<=[MIDNSHP])", perl = TRUE))
    
    current_pos <- read_start
    for (op in ops) {
        op_length <- as.numeric(substr(op, 1, nchar(op) - 1))
        op_type <- substr(op, nchar(op), nchar(op))
        
        if (op_type %in% c("M", "D", "N")) {
            # These operations consume the reference
            if (op_type == "D") {
                # Check if the deletion overlaps with the target range
                deletion_start <- current_pos
                deletion_end <- current_pos + op_length - 1
                if (deletion_start <= target_end && deletion_end >= target_start) {
                    return(TRUE)
                }
            }
            current_pos <- current_pos + op_length
        } else if (op_type %in% c("I", "S", "H")) {
            next
        }
    }
    return(FALSE)
}



findMutationCell <- \(bamFile, bamIndex, 
                      start, end, chr, 
                      reference, alternative,
                      annotate=FALSE,
                      cellBarcode=NULL,
                      cellType=NULL) {
    region <- GRanges(chr, IRanges(start=start, end=end))
    param <- ScanBamParam(which = region, 
                          what = c("seq"),
                          tag = "CB",
                          flag=scanBamFlag(isSecondaryAlignment = FALSE))
    ga <- readGAlignments(bamFile, bamIndex, param=param)
    
    if (nchar(reference)==nchar(alternative)) {
        mcols(ga)$alt <- stackStringsFromGAlignments(ga,region)
        ga <- as.data.table(ga)
        ga <- ga[alt != "."]
        ga[,ref:=reference]
        ga[,nAlt:=sum(alt==alternative), by = CB]
    } else if (nchar(alternative) > nchar(reference)) {
        
        mcols(ga)$ref <- stackStringsFromGAlignments(ga,region)
        ga <- as.data.table(ga)
        ga <- ga[ref != "."]
        target_pos <- start
        ga[,ifAlt:=mapply(.ifInsert, cigar, start, end, target_pos)]
        ga[,nAlt:=sum(ifAlt),by=CB]

    } else if (nchar(alternative) < nchar(reference)) {
        
        mcols(ga)$alt <- stackStringsFromGAlignments(ga,region)
        ga <- as.data.table(ga)
        ga <- ga[alt != "."]
        ga[,ref:=reference]
        target_start <- start
        target_end <- end
        ga[, ifAlt := mapply(.ifDelete, cigar, start, end, 
                             target_start, target_end)]
        ga[,nAlt:=sum(ifAlt),by=CB]
        
    }
    
    ga[,read_depth:=.N, by = CB]
    ga[,alt:=alternative]
    ga <- unique(ga[, .(seqnames, start, end, ref, alt, 
                        read_depth, CB, nAlt)])
    ga[,VAF:=nAlt/read_depth, by = CB]
    ga <- unique(ga[, .(read_depth, CB, nAlt, VAF)])
    if (annotate) {
        idx <- match(ga$CB, cellBarcode)
        ga$cellType <- as.character(cellType[idx])
        ga$cellType[is.na(ga$cellType)] <- "discarded"
    }
    #res <- unique(ga[,.(CB, read_depth, ref, alt, nAlt, VAF, cellType)])
    ga[, `:=`(seqnames = chr, start = start, end = end)]
    ga
}


