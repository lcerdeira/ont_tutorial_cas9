suppressWarnings(
  suppressMessages(
    suppressPackageStartupMessages(
      {
        library(yaml)
        library(parallel)
        library(ShortRead)
        library(pbmcapply)
        library(dplyr)
        library(GenomicRanges)
      }       
    )
  )
)        


config <- yaml.load_file("config.yaml")
bed_src <- config$target_regions
study <- config$study_name
reference <- config$reference_genome
gstride <- as.integer(config$gstride)
target_proximity <- as.integer(config$target_proximity)
offtarget_level <- as.integer(config$offtarget_level)
max_threads <- as.integer(config$threads)

r_results <- file.path("Analysis","R")
on_target <- file.path("Analysis","OnTarget")
off_target <- file.path("Analysis","OffTarget")
dir.create(r_results, showWarnings = FALSE, recursive=TRUE)
dir.create(on_target, showWarnings = FALSE, recursive=TRUE)

##### no more config beyond this point - please #####

referenceFile <- file.path("ReferenceData", basename(reference))

##### place for the common accessory methods ####

loadReferenceGenome <- function() {
  # derive a referenceGenome object from the named fasta elements in the provided fasta reference resource
  referenceGenomeSequence <<- readDNAStringSet(referenceFile)
  referenceGenome <<- data.frame(t(
    data.frame(strsplit(names(referenceGenomeSequence), " "), stringsAsFactors = FALSE)), 
    stringsAsFactors = FALSE)
  referenceGenome$sid <<- seq(nrow(referenceGenome))
  # clip out the mitochondrial genome ...
  if ("MT" %in% referenceGenome[,1]) { # extra lookup to handle cases when reference is a single entry ...
    referenceGenome <<- referenceGenome[-which(referenceGenome[,1]=="MT"),]
  }
}

getStringSetId <- function(chrId) {
  if (!(exists("referenceGenome") & exists("referenceGenomeSequence"))) { loadReferenceGenome() }
  return(referenceGenome[match(as.character(chrId), as.character(referenceGenome[,1])),"sid"])
}

phredmean <- function(l) {
  -10 * log10(mean(10^(l/-10)))
}

qualToMeanQ <- function(qstr) {
  baseq <- as.numeric(charToRaw(qstr))-33
  meanerror <- mean(10^(baseq/-10))
  -10*log10(meanerror)
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
    meanQ <- unlist(lapply(qual.data$qual, qualToMeanQ))
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
cat(paste0("loading reference genome data", "\n"))
loadReferenceGenome()
# for the tutorial we are stripping out the non-primary chromosome chunks ... - in current human genome anything with a dot in name
# the loadReferenceGenome method is removing the MT genome
if (length(grep("\\.", referenceGenome[,1])) > 0) {
  referenceGenome <- referenceGenome[-grep("\\.", referenceGenome[,1]),]
}

####################
# identify the project BAM files
####################
cat(paste0("finding BAM files and loading unmapped reads", "\n"))
mappedBam <- file.path("Analysis","Minimap2", paste0(study, ".bam"))
unmappedBam = file.path("Analysis","Minimap2", paste0(study, ".unmapped.quals"))
# and parse the unmapped content
unmappedReads <- harvestUnmapped(unmappedBam)


####################
# load the bed file
# define what is ontarget, targetproximal and background
####################
cat(paste0("loading BED files and scoring coordinates", "\n"))
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
cat(paste0("loading BAM files", "\n"))
seqlengths(gr) <- width(referenceGenomeSequence[getStringSetId(names(seqlengths(gr)))])
cat(paste0("primary mappings", "\n"))
params=ScanBamParam(which=gr,
                    what=c("mapq", "qual", "qname"), 
                    tag=c("NM"),
                    flag=scanBamFlag(isSupplementaryAlignment=FALSE, isSecondaryAlignment=FALSE))  
wga <- readGAlignments(file=mappedBam, param=params)
# cat(paste0("secondary mappings", "\n"))
# wgaSecondary <- readGAlignments(file=mappedBam, 
#                                 param=ScanBamParam(which=gr,
#                                                    flag=scanBamFlag(isSupplementaryAlignment=FALSE, isSecondaryAlignment=TRUE)))
# cat(paste0("supplementary mappings", "\n"))
# wgaSupplementary <- readGAlignments(file=mappedBam, 
#                                     param=ScanBamParam(which=gr, 
#                                                        flag=scanBamFlag(isSupplementaryAlignment=TRUE, isSecondaryAlignment=FALSE)))

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
cat(paste0("hunting off target mappings", "\n"))
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
cat(paste0("preparing mapping characteristics", "\n"))
backgroundUniverse <- GenomicRanges::reduce(backgroundR)
offtargetUniverse <- GenomicRanges::reduce(offR)
#ontargetUniverse <- GenomicRanges::reduce(br)
###### curious problem here; if the target regions are overlapping (e.g. on different strands ...) then the reduce
# will lead to a size collision in the next names(...) steps - we will skip this reduce step on the ontarget to keep
# it intact and as defined
ontargetUniverse <- ontargetUniverse <- GRanges(br)
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
  
  #startF <- which(start(c1)>=start(nd) & strand(c1)=="+")
  #startR <- which(end(c1)<=end(nd) & strand(c1)=="-")
  #endF <- which(end(c1)<=endF(nd) & strand(c1)=="+")
  #endR <- which(start(c1)>=start(nd) & strand(c1)=="-")
  #starts <- unique(append(startF, startR))
  #stops <-unique(append(endF, endR))  
  
  starts <- which(start(c1)>=start(nd))
  c2 <- c1[starts]
  
  seqlevels(c1) <- unique(as.character(seqnames(c1)))
  seqnames(c1) <- factor(seqnames(c1))
  cov <- GenomicAlignments::coverage(c1, shift=-start(nd), width=width(nd))
  depths <- rep(as.integer(unlist(runValue(cov))), as.integer(unlist(runLength(cov))))
  
  qdata <- quantile(depths)
  
  lf <- letterFrequency(subseq(referenceGenomeSequence[[getStringSetId(seqnames(nd))]], start(nd), end(nd)), c("G", "C", "N"))
  
  rstart <- length(c2)                                       # this is equiv. to earlier readStarts
  basesstart <- sum(qwidth(c2))                              # earlier basesReadsStarted
  meanreadlen <- mean(qwidth(c1))                            # see previous *readLen* - now mean for *any* overlapping read
  startreadlen <- mean(qwidth(c2))                           # what is the mean sequence length for reads that start in segment?
  strandp <- length(which(strand(c1)=="+"))                  # within the interval; not linked to readStarts
  strandn <- length(which(strand(c1)=="-"))
  gccount <- lf[1]+lf[2]
  ncount <- lf[3]
  #mapq <- -10 * log10(mean(10^(-mcols(c2)$mapq/10)))         # mapq for reads starting in segment
  #map0 <- -10 * log10(mean(10^(-mcols(c1)$mapq/10)))         # mapq for reads overlapping segment
  #readq <- -10 * log10(mean(10^(-mean(alphabetScore(mcols(c2)$qual) / width(mcols(c2)$qual))/10))) # per read q score for reads starting in segment 
  #read0 <- -10 * log10(mean(10^(-mean(alphabetScore(mcols(c1)$qual) / width(mcols(c1)$qual))/10))) # per read q score for reads overlapping segment
  
  # discussion with Olle on how Qvalues are best summarised and prepared ... ShortRead alphabetScore is not ideal
  # while alphabetScore(mcols(c2)$qual) / width(mcols(c2)$qual)) gives a number this is not scaled appropriately
  mapq <- phredmean(mcols(c2)$mapq)
  map0 <- phredmean(mcols(c1)$mapq)
  readq <- phredmean(unlist(lapply(as.character(mcols(c2)$qual), qualToMeanQ)))
  read0 <- phredmean(unlist(lapply(as.character(mcols(c1)$qual), qualToMeanQ)))
  
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
  return(c(rstart, basesstart, meanreadlen, startreadlen, strandp, strandn, gccount, ncount, mapq, map0, readq, read0, d005, d025, d050, d075, d095, dmean, nm, cigardel, cigarins, cigarmapped))
}


bamMineUniverse <- function(universe) {
  startStrand <-matrix(unlist(pbmclapply(seq_along(seqnames(universe)), getStartStrand, gdata=universe, mc.cores=min(detectCores()-1, max_threads))), ncol=22, byrow=TRUE)
  universe$rstart <- startStrand[,1]
  universe$basesstart <- startStrand[,2]
  universe$meanreadlen <- startStrand[,3]
  universe$startreadlen <- startStrand[,4]
  universe$strandp <- startStrand[,5]
  universe$strandn <- startStrand[,6]
  universe$gccount <- startStrand[,7]
  universe$ncount <- startStrand[,8]
  universe$mapq <- startStrand[,9]
  universe$map0 <- startStrand[,10]
  universe$readq <- startStrand[,11]
  universe$read0 <- startStrand[,12]
  universe$d005 <- startStrand[,13]
  universe$d025 <- startStrand[,14]
  universe$d050 <- startStrand[,15]
  universe$d075 <- startStrand[,16]
  universe$d095 <- startStrand[,17]
  universe$dmean <- startStrand[,18]
  universe$nm <- startStrand[,19]
  universe$cigardel <- startStrand[,20]
  universe$cigarins <- startStrand[,21]
  universe$cigarmapped <- startStrand[,22]
  return(universe)
}


####################
# and apply the function to the mapping universe types ...
####################
cat(paste0("background", "\n"))
backgroundUniverse <- bamMineUniverse(backgroundUniverse)
cat(paste0("offtarget", "\n"))
offtargetUniverse <- bamMineUniverse(offtargetUniverse)
cat(paste0("ontarget", "\n"))
ontargetUniverse <- bamMineUniverse(ontargetUniverse)
cat(paste0("target proximal", "\n"))
targetproximalUniverse <- bamMineUniverse(targetproximalUniverse)


####################
# and save the data in a file for import into the report ...
####################
mappingResultsFile <- file.path(r_results, paste0(study, "_mapping_results", ".Rdata"))
save(br, wga.cov, backgroundUniverse, offtargetUniverse, ontargetUniverse, targetproximalUniverse, file=mappingResultsFile)

####################
# create the finer resolution aggregate data for the discrete target regions ...
####################
aggregateDepthInfo <- function(x, xr, ontarget=FALSE, geneId=NA) {
  targetRegion <- xr[x]
  proximalRegion <- targetRegion
  if (ontarget) {
    start(proximalRegion) <- max(start(proximalRegion) - target_proximity, 1)
    end(proximalRegion) <- end(proximalRegion) + target_proximity
  }
  if (is.na(geneId)) { 
    geneId <- names(targetRegion) 
  }
  if (geneId != "OffTarget") {
    cat(paste0("geneId:", geneId, "\n"))
  }
  
  c0 <- wga[seqnames(wga)==as.character(seqnames(targetRegion))]
  c1 <- c0[subjectHits(findOverlaps(proximalRegion, granges(c0)))]
  c2 <- c0[subjectHits(findOverlaps(targetRegion, granges(c0)))]
  seqlevels(c1) <- unique(as.character(seqnames(c1)))
  seqnames(c1) <- factor(seqnames(c1))
  cov <- GenomicAlignments::coverage(c1, shift=-start(proximalRegion), width=width(proximalRegion))

  seqlen <- width(proximalRegion)
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
  fcov <- GenomicAlignments::coverage(c1[which(strand(c1)=="+")], shift=-start(proximalRegion), width=width(proximalRegion))
  ba$fwd_cov <- mcols(binnedAverage(bins, fcov, "fwd_cov"))$fwd_cov
  
  # write the target sequences to file ...
  if (ontarget) {
    write(mcols(c2)$qname, file.path(on_target, paste0(study, ".", geneId, ".mappedreads")), append=TRUE)
  } else {
    write(mcols(c2)$qname, file.path(off_target, paste0(study, ".", geneId, ".mappedreads")), append=TRUE)
  }
  return(as.data.frame(ba))
}
  


cat(paste0("scoring data per target region", "\n"))
suppressWarnings({
  aggregatedCov <- bind_rows(pbmclapply(seq_along(ontargetUniverse), aggregateDepthInfo, xr=ontargetUniverse, ontarget=TRUE, mc.cores=min(detectCores()-1, max_threads)), .id = "column_label")
  aggregatedCovFile <- file.path(r_results, paste0(study, "_aggregated_coverage", ".Rdata"))
  aggregatedGR <- GenomicRanges::makeGRangesFromDataFrame(aggregatedCov[,-1], keep.extra.columns = TRUE)
  save(aggregatedGR, file=aggregatedCovFile)
})


cat(paste0("parsing off-target/background coverage - please be patient ...", "\n"))
suppressWarnings({
  offtCov <- pbmclapply(seq_along(offtargetUniverse), aggregateDepthInfo, xr=offtargetUniverse, geneId="OffTarget", mc.cores=min(detectCores()-1,max_threads))
  aggregatedOff <- bind_rows(offtCov, .id = "column_label")
  aggregatedOffFile <- file.path(r_results, paste0(study, "_aggregated_offt_coverage", ".Rdata"))
  save(aggregatedOff, file=aggregatedOffFile)
})


