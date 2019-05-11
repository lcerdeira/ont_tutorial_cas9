library(yaml)
library(parallel)
library(ShortRead)
library(pbmcapply)
library(dplyr)
library(GenomicRanges)

config <- yaml.load_file("config.yaml")
bed_src <- config$target_regions
study <- config$study_name
reference <- config$reference_genome
gstride <- as.integer(config$gstride)
target_proximity <- as.integer(config$target_proximity)
offtarget_level <- as.integer(config$offtarget_level)

r_results <- file.path("Analysis","R")
on_target <- file.path("Analysis","OnTarget")
dir.create(r_results, showWarnings = FALSE, recursive=TRUE)
dir.create(on_target, showWarnings = FALSE, recursive=TRUE)

##### no more config beyond this point - please #####

##### place for the common accessory methods ####

loadReferenceGenome <- function() {
  # derive a referenceGenome object from the named fasta elements in the provided fasta reference resource
  referenceGenomeSequence <<- readDNAStringSet(referenceGenome.file)
  referenceGenome <<- data.frame(t(
    data.frame(strsplit(names(referenceGenomeSequence), " "), stringsAsFactors = FALSE)), 
    stringsAsFactors = FALSE)
  referenceGenome$sid <<- seq(nrow(referenceGenome))
  # clip out the mitochondrial genome ...
  referenceGenome <<- referenceGenome[-which(referenceGenome[,1]=="MT"),]
}

getStringSetId <- function(chrId) {
  if (!(exists("referenceGenome") & exists("referenceGenomeSequence"))) { loadReferenceGenome() }
  return(referenceGenome[match(as.character(chrId), as.character(referenceGenome[,1])),"sid"])
}

harvestUnmapped <- function(qualfilelocation, chunk.size=100000, force=FALSE) {
  # time samtools view -@ 5 -f 0x4 -O BAM clive_ngmlr.samtools.bam > unmapped.bam
  # samtools view -@ 5 -O sam unmapped.bam | awk '{print $11}' > unmapped.quals
  
  # samtools view -@ 4 -O sam -h clive.minimap2.bam | samtools view -@ 4 -F 0x4 -O BAM -U clive.minimap2.unmapped.bam -T ../../ReferenceData/Homo_sapiens.GRCh37.75.dna.chromosome.ALL.fasta | samtools sort -@ 4 -T ../../ReferenceData/Homo_sapiens.GRCh37.75.dna.chromosome.ALL.fasta > clive.minimap2.samtools.bam
  
  chromosomeFile <- file.path(r_results, paste(sub("\\.[^.]*$", "", basename(qualfilelocation)), "rcounts", "Rdata",sep="."))
  if (file.exists(chromosomeFile) & !force) {
    unmapped.content <- readRDS(file=chromosomeFile)
    return(unmapped.content)
  }
  offset <- 0
  unmapped.content <- data.frame(width=integer(),
                                 quality=numeric(),
                                 stringsAsFactors=FALSE) 
  repeat {
    cat(paste("iter", offset, "\n"))
    qual.data <- data.table::fread(file=qualfilelocation, nrows=chunk.size, skip=offset, sep="\n", col.names=c("qual"))
    width <- nchar(qual.data$qual)
    meanQ <- alphabetScore(FastqQuality(qual.data$qual)) / width
    unmapped.content <- rbind(unmapped.content, data.frame(width=width, quality=meanQ))
    if (nrow(qual.data) < chunk.size) {
      break
    }
    offset <- offset + chunk.size
  }
  saveRDS(unmapped.content, file=chromosomeFile)
  return(unmapped.content)
}


##### end of accessory methods #####

####################
# load the reference genome
####################
referenceGenome.file <- file.path("Analysis", "ReferenceGenome", gsub("(\\.zip)|(\\.bz2)|(\\.gz)", "", basename(reference)))
loadReferenceGenome()
# for the tutorial we are stripping out the non-primary chromosome chunks ... - in current human genome anything with a dot in name
# the loadReferenceGenome method is removing the MT genome
referenceGenome <- referenceGenome[-grep("\\.", referenceGenome[,1]),]


####################
# identify the project BAM files
####################
mappedBam <- file.path("Analysis","Minimap2", paste0(study, ".bam"))
unmappedBam = file.path("Analysis","Minimap2", paste0(study, ".unmapped.quals"))
# and parse the unmapped content
unmappedReads <- harvestUnmapped(unmappedBam)


####################
# load the bed file
# define what is ontarget, targetproximal and background
####################
bed <- data.table::fread(file=bed_src)
# create a genomic ranges object for the bed file elements
br <- GRanges(seqnames=unlist(bed[,1]), IRanges(start=unlist(bed[,2]), end=unlist(bed[,3])))
# define flanking regions for the target-proximal analysis
fr <- union(flank(br, width=target_proximity, start=TRUE), flank(br, width=target_proximity, start=FALSE))
# create a genomic ranges object for the chromosomes
gr <- GRanges(seqnames=referenceGenome[,1], IRanges(start=1, end=nchar(referenceGenomeSequence)[referenceGenome[,5]]))
# and define an off-target ranges ..
or <- setdiff(disjoin(c(gr, GenomicRanges::reduce(union(br, fr)))), GenomicRanges::reduce(union(br, fr)))
# name the genes associated with the BED annotations for subsequent display ...
names(br) <- unlist(bed[,4])


####################
# load the BAM file into memory; parse depths of coverage
# STRIDE from YAML; 
####################
seqlengths(gr) <- width(referenceGenomeSequence[getStringSetId(names(seqlengths(gr)))])
params=ScanBamParam(which=gr,
                    what=c("mapq", "qual", "qname"), 
                    tag=c("NM"),
                    flag=scanBamFlag(isSupplementaryAlignment=FALSE, isSecondaryAlignment=FALSE))  
wga <- readGAlignments(file=mappedBam, param=params)

wgaSecondary <- readGAlignments(file=mappedBam, 
                                param=ScanBamParam(which=gr, 
                                                   flag=scanBamFlag(isSupplementaryAlignment=FALSE, isSecondaryAlignment=TRUE)))
wgaSupplementary <- readGAlignments(file=mappedBam, 
                                    param=ScanBamParam(which=gr, 
                                                       flag=scanBamFlag(isSupplementaryAlignment=TRUE, isSecondaryAlignment=FALSE)))

max.read.len <- max(width(wga))
wgac <- coverage(wga)
bins <- tileGenome(seqlengths(gr), tilewidth=gstride, cut.last.tile.in.chrom=TRUE)
# handle any hanger-on chromosome ids
wgac <- wgac[which(names(wgac) %in% seqlevels(bins))]
wgaba.all <- binnedAverage(bins, wgac, "binned_cov")


####################
# look for peaks of potential off-target; threshold from YAML
####################
# strip out on target and flanking regions ...
wgaba <- wgaba.all[-queryHits(findOverlaps(wgaba.all, br))]
wgaba <- wgaba[-queryHits(findOverlaps(wgaba, fr))]
wga.cov <- mean(wgaba$binned_cov) # calculate mean background coverage

offR <- wgaba[which(wgaba$binned_cov > (offtarget_level*wga.cov))]
backgroundR <- wgaba[-queryHits(findOverlaps(wgaba, offR))]
#backgroundRDepth <- backgroundR[which(backgroundR$binned_cov > 0)]



####################
# Reduce the background dimensions - not aiming to plot chromosomal coverages ...
# this is available in other tutorials
# selecting simpler names for the objects; they will be reused ...
####################
backgroundUniverse <- GenomicRanges::reduce(backgroundR)
offtargetUniverse <- GenomicRanges::reduce(offR)
ontargetUniverse <- GenomicRanges::reduce(br)
targetproximalUniverse <- GenomicRanges::reduce(fr)

names(ontargetUniverse) <- names(br)


####################
# big ugly BAM parsing function; dive across a GRange
####################
getStartStrand <- function(x, gdata) {
  nd <- gdata[x]
  # how to perform a fast reduction on data size - filter by chromosome
  c0 <- wga[seqnames(wga)==as.character(seqnames(nd))]
  c1 <- c0[subjectHits(findOverlaps(nd, granges(c0)))]
  c2 <- c1[start(c1)>=start(nd)]
  c3 <- c1[end(c1)<=end(nd)]
  
  sc0 <- wgaSecondary[seqnames(wgaSecondary)==as.character(seqnames(nd))]
  sc1 <- sc0[subjectHits(findOverlaps(nd, granges(sc0)))]
  sc2 <- sc1[start(sc1)>=start(nd)]
  
  sp0 <- wgaSupplementary[seqnames(wgaSupplementary)==as.character(seqnames(nd))]
  sp1 <- sp0[subjectHits(findOverlaps(nd, granges(sp0)))]
  sp2 <- sp1[start(sp1)>=start(nd)]
  
  seqlevels(c1) <- unique(as.character(seqnames(c1)))
  seqnames(c1) <- factor(seqnames(c1))
  cov <- GenomicAlignments::coverage(c1, shift=-start(nd), width=width(nd))
  depths <- rep(as.integer(unlist(runValue(cov))), as.integer(unlist(runLength(cov))))
  qdata <- quantile(depths)
  
  lf <- letterFrequency(subseq(referenceGenomeSequence[[getStringSetId(seqnames(nd))]], start(nd), end(nd)), c("G", "C", "N"))
  
  rstart <- length(c2)                                        # this is equiv. to earlier readStarts
  basesstart <- sum(qwidth(c2))                              # earlier basesReadsStarted
  meanreadlen <- mean(qwidth(c1))                            # see previous *readLen* - now mean for *any* overlapping read
  startreadlen <- mean(qwidth(c2))                           # what is the mean sequence length for reads that start in segment?
  endreadlen <- mean(qwidth(c3))                             # what is the mean sequence length for reads that end in the segment?
  strandp <- length(which(strand(c1)=="+"))                  # within the interval; not linked to readStarts
  strandn <- length(which(strand(c1)=="-"))
  gccount <- lf[1]+lf[2]
  ncount <- lf[3]
  mapq <- mean(mcols(c2)$mapq)                               # mapq for reads starting in segment
  map0 <- mean(mcols(c1)$mapq)                               # mapq for reads overlapping segment
  readq <- mean(alphabetScore(mcols(c2)$qual) / width(mcols(c2)$qual)) # per read q score for reads starting in segment 
  read0 <- mean(alphabetScore(mcols(c1)$qual) / width(mcols(c1)$qual)) # per read q score for reads overlapping segment
  d005 <- qdata[1]
  d025 <- qdata[2]
  d050 <- qdata[3]
  d075 <- qdata[4]
  d095 <- qdata[5]
  dmean <- mean(depths)
  nm <- sum(mcols(c2)$NM)                                    # this is the #NM mismatch count; reads starting in segment
  cigardel <- sum(sum(width(cigarRangesAlongPairwiseSpace(cigar(c2), ops=c("D")))))
  cigarins <- sum(sum(width(cigarRangesAlongPairwiseSpace(cigar(c2), ops=c("I")))))
  cigarmapped <- sum(cigarWidthAlongQuerySpace(cigar(c2), after.soft.clipping=TRUE))
  isseccount <- length(sc2)
  issecbases <- sum(qwidth(sc2))
  issuppcount <- length(sp2)
  issuppbases <- sum(qwidth(sp2))
  return(c(rstart, basesstart, meanreadlen, startreadlen, endreadlen, strandp, strandn, gccount, ncount, mapq, map0, readq, read0, d005, d025, d050, d075, d095, dmean, nm, cigardel, cigarins, cigarmapped, isseccount, issecbases, issuppcount, issuppbases))
}


bamMineUniverse <- function(universe) {
  startStrand <-matrix(unlist(pbmclapply(seq_along(seqnames(universe)), getStartStrand, gdata=universe, mc.cores=min(detectCores()-1,8))), ncol=27, byrow=TRUE)
  universe$rstart <- startStrand[,1]
  universe$basesstart <- startStrand[,2]
  universe$meanreadlen <- startStrand[,3]
  universe$startreadlen <- startStrand[,4]
  universe$endreadlen <- startStrand[,5]
  universe$strandp <- startStrand[,6]
  universe$strandn <- startStrand[,7]
  universe$gccount <- startStrand[,8]
  universe$ncount <- startStrand[,9]
  universe$mapq <- startStrand[,10]
  universe$map0 <- startStrand[,11]
  universe$readq <- startStrand[,12]
  universe$read0 <- startStrand[,13]
  universe$d005 <- startStrand[,14]
  universe$d025 <- startStrand[,15]
  universe$d050 <- startStrand[,16]
  universe$d075 <- startStrand[,17]
  universe$d095 <- startStrand[,18]
  universe$dmean <- startStrand[,19]
  universe$nm <- startStrand[,20]
  universe$cigardel <- startStrand[,21]
  universe$cigarins <- startStrand[,22]
  universe$cigarmapped <- startStrand[,23]
  universe$isseccount <- startStrand[,24]
  universe$issecbases <- startStrand[,25]
  universe$issuppcount <- startStrand[,26]
  universe$issuppbases <- startStrand[,27]
  return(universe)
}


####################
# and apply the function to the mapping universe types ...
####################
backgroundUniverse <- bamMineUniverse(backgroundUniverse)
offtargetUniverse <- bamMineUniverse(offtargetUniverse)
ontargetUniverse <- bamMineUniverse(ontargetUniverse)
targetproximalUniverse <- bamMineUniverse(targetproximalUniverse)


####################
# and save the data in a file for import into the report ...
####################
mappingResultsFile <- file.path(r_results, paste0(study, "_mapping_results", ".Rdata"))
save(br, backgroundUniverse, offtargetUniverse, ontargetUniverse, targetproximalUniverse, file=mappingResultsFile)

####################
# create the finer resolution aggregate data for the discrete target regions ...
####################
aggregateDepthInfo <- function(x, xr, geneId=NA) {
  targetRegion <- xr[x]
  if (is.na(geneId)) { 
    geneId <- names(targetRegion) 
  }
  cat(paste0("geneId:", geneId, "\n"))
  
  c0 <- wga[seqnames(wga)==as.character(seqnames(targetRegion))]
  c1 <- c0[subjectHits(findOverlaps(targetRegion, granges(c0)))]
  seqlevels(c1) <- unique(as.character(seqnames(c1)))
  seqnames(c1) <- factor(seqnames(c1))
  cov <- GenomicAlignments::coverage(c1, shift=-start(targetRegion), width=width(targetRegion))

  seqlen <- width(targetRegion)
  names(seqlen)<-as.character(seqnames(targetRegion))
  bins <- GenomicRanges::tileGenome(seqlengths=seqlen, tilewidth=10, cut.last.tile.in.chrom=TRUE)
  ba <- binnedAverage(bins, cov, "binned_cov")
  slen <- width(referenceGenomeSequence[getStringSetId(unique(seqnames(ba)))])
  names(slen) <- unique(seqnames(ba))
  seqlengths(ba) <- slen

  ba$pos <- start(ba)
  end(ba) <- end(ba)+start(targetRegion)
  start(ba) <- start(ba)+start(targetRegion)
  ba$gene <- geneId
  
  # calculate the coverage of reads on fwd strand for directionality plotting
  fcov <- GenomicAlignments::coverage(c1[which(strand(c1)=="+")], shift=-start(targetRegion), width=width(targetRegion))
  ba$fwd_cov <- mcols(binnedAverage(bins, fcov, "fwd_cov"))$fwd_cov
  
  # write the target sequences to file ...
  write(mcols(c1)$qname, file.path(on_target, paste0(geneId, ".mappedreads")), append=TRUE)
  
  return(as.data.frame(ba))
}
  



aggregatedCov <- bind_rows(pbmclapply(seq_along(ontargetUniverse), aggregateDepthInfo, xr=ontargetUniverse, mc.cores=min(detectCores()-1,8)), .id = "column_label")
aggregatedCovFile <- file.path(r_results, paste0(study, "_aggregated_coverage", ".Rdata"))

aggregatedGR <- GenomicRanges::makeGRangesFromDataFrame(aggregatedCov[,-1], keep.extra.columns = TRUE)

save(aggregatedGR, file=aggregatedCovFile)


offtCov <- pbmclapply(seq_along(offtargetUniverse), aggregateDepthInfo, xr=offtargetUniverse, geneId="OffTarget", mc.cores=min(detectCores()-1,8))
