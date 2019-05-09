library(session)
source("BamParser.R")

args = commandArgs(trailingOnly=TRUE)
referenceGenome.file <- args[1]
mappedBam <- args[2]
unmappedBam <- args[3]

r_results <- file.path("Analysis","R")
dir.create(r_results, showWarnings = FALSE, recursive=TRUE)

cat(paste("ReferenceGenome :", referenceGenome.file, "\n"))
cat(paste("mappedBam       :", mappedBam, "\n"))
cat(paste("unmappedBam     :", unmappedBam, "\n"))

loadReferenceGenome()
# for the tutorial we are stripping out the non-primary chromosome chunks ...
referenceGenome <- referenceGenome[-grep("\\.", referenceGenome[,1]),]

# parseBamFile(mappedBam)
# harvestUnmapped(unmappedBam)




bed <- data.table::fread(file=file.path("RawData", "enrichment_targets.bed"))
# create a genomic ranges object for the bed file elements
br <- GRanges(seqnames=unlist(bed[,1]), IRanges(start=unlist(bed[,2]), end=unlist(bed[,3])))
# define flanking regions for the target-proximal analysis
flanking.size <- 25000
fr <- union(flank(br, width=flanking.size, start=TRUE), flank(br, width=flanking.size, start=FALSE))
# create a genomic ranges object for the chromosomes
gr <- GRanges(seqnames=referenceGenome[,1], IRanges(start=1, end=nchar(referenceGenomeSequence)[referenceGenome[,5]]))
# and define an off-target ranges ..
or <- setdiff(disjoin(c(gr, reduce(union(br, fr)))), reduce(union(br, fr)))
# name the genes associated with the BED annotations for subsequent display ...
names(br) <- unlist(bed[,4])

## calculate mapping characteristics for or (off target), fr (flanking regions) and br (on target regions)
extractMappingInfo <- function(x, xr) {
  chrId <- as.character(seqnames(xr[x]))
  start <- start(ranges(xr[x]))
  window <- width(ranges(xr[x]))
  sid <- getStringSetId(chrId)
  harvestBam(start, sid, chrId, window, mappedBam) 
}
orRes <- fixBamFileColumns(as.data.frame(t(as.data.frame(lapply(seq_along(or), extractMappingInfo, xr=or))), stringsAsFactors=FALSE, row.names=seq_along(or)))

frRes <- fixBamFileColumns(as.data.frame(t(as.data.frame(lapply(seq_along(fr), extractMappingInfo, xr=fr))), stringsAsFactors=FALSE, row.names=seq_along(fr)))

brRes <- fixBamFileColumns(as.data.frame(t(as.data.frame(lapply(seq_along(br), extractMappingInfo, xr=br))), stringsAsFactors=FALSE, row.names=seq_along(br)))

