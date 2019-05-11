
library(Rsamtools)
library(parallel)
library(GenomicAlignments)
library(dplyr)
library(gtools)
library(pbmcapply)
library(ShortRead)



# prepare an accessory method for mapping named chromosomes to their pointers in the reference fasta
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


quickHarvest <- function(x, y, chrId, bamfilelocation, supplementary, secondary) {
  nsupp <- 0
  nsecd <- 0
  par=ScanBamParam(which=GRanges(seqnames = chrId, ranges = IRanges(start = x, end = y)), 
                      what=c("strand", "pos", "qwidth", "mapq", "cigar"), 
                      flag=scanBamFlag(isSupplementaryAlignment=supplementary, isSecondaryAlignment=secondary), 
                      tag=c("NM"))
  cigar.df <- as.data.frame(scanBam(BamFile(bamfilelocation), param=par)[[1]])
  ons <- which(cigar.df$pos >= x)
  cigar.m <- cigarWidthAlongQuerySpace(cigar.df[ons, "cigar"], after.soft.clipping=TRUE)
  cigar.d <- width(cigarRangesAlongPairwiseSpace(cigar.df[ons, "cigar"], ops=c("D")))
  cigar.i <- width(cigarRangesAlongPairwiseSpace(cigar.df[ons, "cigar"], ops=c("I")))
  cigar.d.bases <- unlist(lapply(cigar.d, sum))
  cigar.i.bases <- unlist(lapply(cigar.i, sum))
  meanCov <- 0
  if (nrow(cigar.df)>0) {
    cov <- GAlignments(seqnames=rep("map", nrow(cigar.df)), pos=cigar.df$pos, strand=cigar.df$strand, cigar=as.character(cigar.df$cigar))
    covTable <- as.data.frame(coverage(cov, shift=-x, width=(y-x+1)))
    meanCov <- mean(covTable$value)
  }
  if (supplementary) { nsupp <- length(ons) }
  if (secondary) { nsecd <- length(ons) }
  return(c(xsupp=nsupp,
           xsecd=nsecd,
           xcov=round(meanCov,digits=2), 
           xmapped=sum(cigar.m), 
           xmapq=round(mean(cigar.df[ons,"mapq"]), digits=2),
           xwidth=round(mean(cigar.df[ons,"qwidth"]), digits=2),
           xmismatch=sum(cigar.df[ons,"NM"]), 
           xins=sum(cigar.i.bases), 
           xdel=sum(cigar.d.bases)))
}


harvestBam <- function(x, dnaStringSetId, chrId, window.size, bamfilelocation) {
  y = x+window.size
  x <- x + 1
  if (y > nchar(referenceGenomeSequence[[dnaStringSetId]])) {
    y = nchar(referenceGenomeSequence[[dnaStringSetId]])
  }

  params=ScanBamParam(which=GRanges(seqnames = chrId, ranges = IRanges(start = x, end = y)), 
                      what=c("flag", "strand", "pos", "qwidth", "mapq", "cigar", "qual"), 
                      flag=scanBamFlag(isSupplementaryAlignment=FALSE, isSecondaryAlignment=FALSE), 
                      tag=c("NM"))
  SeqCigar <- as.data.frame(scanBam(BamFile(bamfilelocation), param=params)[[1]])
  
  # to ensure for single counting of reads - select for reads that start within the defined window
  # rather than all reads overlapping the defined window
  onStart <- which(SeqCigar$pos >= x)
  
  # enumerate the base space which has been mapped from FASTQ
  cigarMapped <- cigarWidthAlongQuerySpace(SeqCigar[onStart, "cigar"], after.soft.clipping=TRUE)

  # count CIGAR derived INS and DEL events
  cigarDEL <- width(cigarRangesAlongPairwiseSpace(SeqCigar[onStart, "cigar"], ops=c("D")))
  cigarINS <- width(cigarRangesAlongPairwiseSpace(SeqCigar[onStart, "cigar"], ops=c("I")))
  cigarDELevents <- unlist(lapply(cigarDEL, length))
  cigarINSevents <- unlist(lapply(cigarINS, length))
  cigarDELbases <- unlist(lapply(cigarDEL, sum))
  cigarINSbases <- unlist(lapply(cigarINS, sum))
  
    
  if (!(exists("referenceGenome") & exists("referenceGenomeSequence"))) { loadReferenceGenome() }
  letterFreq <- letterFrequency(subseq(referenceGenomeSequence[[dnaStringSetId]], x, y), c("A", "C", "G", "T", "N"))
  width=sum(letterFreq)
  
  # parse out the mean per-read base calling qvalues
  alphabetScore <- alphabetScore(FastqQuality(SeqCigar[onStart, "qual"])) / SeqCigar[onStart, "qwidth"]
  
  # process the depth of coverage stuff ...
  #tabulatedBases <- tabulate(unlist(apply(SeqCigar[,c("pos", "qwidth")], 1, getBases)))[seq(x, y)]
  #tb <- quantile(tabulatedBases, na.rm=TRUE)
  #names(tb) <- paste("IQR",names(tb))
  tb <- c(0,0,0,0,0)
  meanCov <- 0
  if (nrow(SeqCigar)>0) {
    cov <- GAlignments(seqnames=rep("map", nrow(SeqCigar)), pos=SeqCigar$pos, strand=SeqCigar$strand, cigar=as.character(SeqCigar$cigar))
    covTable <- as.data.frame(coverage(cov, shift=-x, width=window.size))
    tb <- quantile(covTable$value)
    meanCov <- mean(covTable$value)
  }
  names(tb) <- c("iqr.0", "iqr.25", "iqr.50", "iqr.75", "iqr.100")
  
  # and some basic characteristics for the supplemental reads
  # supplFlags <- scanBamFlag(isSupplementaryAlignment=TRUE)
  
  
  res <- c(chrId=chrId,
    startpos=x,
    readStarts=length(onStart),                        # validated as correct with Chr20 dataset (samtools as ground truth)
    plusStrand=length(which(SeqCigar[onStart, "strand"]=="+")),
    basesReadsStarted=sum(SeqCigar[onStart,"qwidth"]), # validated as correct with Chr20 dataset
    width=width,                                       # the width of the window considered
    gccount=as.integer(letterFreq['G'] + letterFreq['C']),
    ncount=as.integer(letterFreq['N']),
    mismatches=sum(SeqCigar[onStart,"NM"]),            # the sum of total mapping edit distance
    cigarMapped=sum(cigarMapped),                      # the sum of read bases that are mapped to reference genome, starting in this window, after softclipping
    cigarInsertionBases=sum(cigarINSbases),            # the sum of inserted bases against the reference genome (pairwise presentation)
    cigarDeletionBases=sum(cigarDELbases),             # the sum of inserted bases against the reference genome (pairwise presentation)
    cigarInsertionEvents=sum(cigarINSevents),
    cigarDeletionEvents=sum(cigarDELevents),
    mapq=round(mean(SeqCigar[onStart,"mapq"]), digits=2),# mean mapping quality for reads starting in this window
    readQ=round(mean(alphabetScore), digits=2),        # this is the per-read mean Q value averaged across onStart
    readLen=round(mean(mean(SeqCigar[onStart, "qwidth"])),digits=2),  # average sequence length for onStart reads
    unlist(tb),                                                 # the depth of coverage summary
    meanCov=round(meanCov, digits=2),
    flagCount=length(unique(SeqCigar$flag)),
    colSums(rbind(quickHarvest(x, y, chrId, bamfilelocation, TRUE, FALSE), quickHarvest(x, y, chrId, bamfilelocation, FALSE, TRUE)), na.rm=TRUE)
  )
  return(res)
}


harvestChromosome <- function(chrId, bamfilelocation, window.size=100000, mc.cores=min(detectCores()-1, 24), force=FALSE) {
  chromosomeFile <- file.path(r_results, paste(sub("\\.[^.]*$", "", basename(bamfilelocation)), paste("chrId",chrId,as.integer(window.size),sep="_"), "Rdata",sep="."))
  if (file.exists(chromosomeFile) & !force) {
     chrData <- readRDS(file=chromosomeFile)
     return(chrData)
  }
  if (!(exists("referenceGenome") & exists("referenceGenomeSequence"))) { loadReferenceGenome() }
  dnaStringSetId <- getStringSetId(chrId)
  coords <- seq(0, length(referenceGenomeSequence[[dnaStringSetId]]), by=window.size)
  cat(paste("harvest - ", chrId, dnaStringSetId, paste(range(coords), collapse=" "), "\n"))
  # harv <- lapply(coords, harvestBam, dnaStringSetId=dnaStringSetId, chrId=chromosomeId)
  # this takes approx 476 seconds for Chr20
  # just pulling in the data from the BAM file takes around 310s - this is uncompress bound?
  # with a new bamfile object for each iteration through lapply time approx 310s
  
  # an attempt with mclapply ...
#mcharv <- lapply(coords, harvestBam, dnaStringSetId=dnaStringSetId, chrId=chrId)
  # using 2 cores this halved the time ... to 154s
  # using all 8 available cores ... to 41s....
mcharv <- pbmclapply(coords, harvestBam, dnaStringSetId=dnaStringSetId, chrId=chrId, window.size=window.size, bamfilelocation=bamfilelocation, mc.cores=mc.cores, mc.preschedule=FALSE, mc.silent=FALSE)
  
  # a little extra output management ... remove any misshapen results ...
  # expectedRows <- length(mcharv[[1]])
  # expectedGood <- which(unlist(lapply(mcharv, length))==expectedRows)

# cl <- makeCluster(5, type="FORK")
# clusterEvalQ(cl, library(Rsamtools))
# clusterEvalQ(cl, library(GenomicAlignments))
# 
# clusterExport(cl, 'bamfile.path')
# clusterExport(cl, 'window.size')
# clusterExport(cl, 's')
# message(date())
# mcharv <- parLapply(cl, coords, harvestBam, dnaStringSetId=dnaStringSetId, chrId=chrId)
# message(date())
# stopCluster(cl)

    chrData <- fixBamFileColumns(as.data.frame(t(as.data.frame(mcharv, stringsAsFactors=FALSE)), 
                           col.names=names(mcharv[[1]]), 
                           row.names=seq_along(length(mcharv)), 
                           stringsAsFactors=FALSE))
    
  saveRDS(chrData, file=chromosomeFile)
  return(chrData)
}


parseBamFile <- function(bamfilelocation, window.size=100000, mc.cores=min(detectCores()-1, 24), force=FALSE) {
  parseBamFileResults <- file.path(r_results, paste(sub("\\.[^.]*$", "", basename(bamfilelocation)), paste("aggregated",as.integer(window.size),sep="_"), "Rdata",sep="."))
  if (file.exists(parseBamFileResults) & !force) {
    parsedBam <- readRDS(file=parseBamFileResults)
    return(parsedBam)
  }
  if (!(exists("referenceGenome") & exists("referenceGenomeSequence"))) { loadReferenceGenome() }
  parsedBam <- bind_rows(lapply(gtools:::mixedsort(referenceGenome[,1]), harvestChromosome, force=force, bamfilelocation=bamfilelocation, window.size=window.size, mc.cores=mc.cores), .id = "column_label")
  saveRDS(parsedBam, file=parseBamFileResults)
  return(parsedBam)
}


# the raw data table loses all data type information - the vector returned from harvestBam is heterogeneous
# and is converted into a data.frame and transposed ... some manual intervention is required 
fixBamFileColumns <- function(parsedBamFile) {
  if (!is.integer(parsedBamFile$startpos)) { parsedBamFile$startpos <- as.integer(parsedBamFile$startpos) }
  if (!is.integer(parsedBamFile$readStarts)) { parsedBamFile$readStarts <- as.integer(parsedBamFile$readStarts) }
  if (!is.integer(parsedBamFile$plusStrand)) { parsedBamFile$plusStrand <- as.integer(parsedBamFile$plusStrand) }
  if (!is.integer(parsedBamFile$basesReadsStarted)) { parsedBamFile$basesReadsStarted <- as.integer(parsedBamFile$basesReadsStarted) }
  if (!is.integer(parsedBamFile$width)) { parsedBamFile$width <- as.integer(parsedBamFile$width) }
  if (!is.integer(parsedBamFile$gccount)) { parsedBamFile$gccount <- as.integer(parsedBamFile$gccount) }
  if (!is.integer(parsedBamFile$ncount)) { parsedBamFile$ncount <- as.integer(parsedBamFile$ncount) }
  if (!is.integer(parsedBamFile$mismatches)) { parsedBamFile$mismatches  <- as.integer(parsedBamFile$mismatches) }
  if (!is.integer(parsedBamFile$cigarMapped)) { parsedBamFile$cigarMapped  <- as.integer(parsedBamFile$cigarMapped) }
  if (!is.integer(parsedBamFile$cigarInsertionBases)) { parsedBamFile$cigarInsertionBases  <- as.integer(parsedBamFile$cigarInsertionBases) }
  if (!is.integer(parsedBamFile$cigarDeletionBases)) { parsedBamFile$cigarDeletionBases  <- as.integer(parsedBamFile$cigarDeletionBases) }
  if (!is.integer(parsedBamFile$cigarInsertionEvents)) { parsedBamFile$cigarInsertionEvents  <- as.integer(parsedBamFile$cigarInsertionEvents) }
  if (!is.integer(parsedBamFile$cigarDeletionEvents)) { parsedBamFile$cigarDeletionEvents  <- as.integer(parsedBamFile$cigarDeletionEvents) }
  if (!is.integer(parsedBamFile$iqr.0)) { parsedBamFile$iqr.0  <- as.integer(parsedBamFile$iqr.0) }
  if (!is.integer(parsedBamFile$iqr.25)) { parsedBamFile$iqr.25  <- as.integer(parsedBamFile$iqr.25) }
  if (!is.integer(parsedBamFile$iqr.50)) { parsedBamFile$iqr.50  <- as.integer(parsedBamFile$iqr.50) }
  if (!is.integer(parsedBamFile$iqr.75)) { parsedBamFile$iqr.75  <- as.integer(parsedBamFile$iqr.75) }
  if (!is.integer(parsedBamFile$iqr.100)) { parsedBamFile$iqr.100  <- as.integer(parsedBamFile$iqr.100) }
  if (!is.integer(parsedBamFile$xsupp)) { parsedBamFile$xsupp  <- as.integer(parsedBamFile$xsupp) }
  if (!is.integer(parsedBamFile$xsecd)) { parsedBamFile$xsecd  <- as.integer(parsedBamFile$xsecd) }
  if (!is.integer(parsedBamFile$xmapped)) { parsedBamFile$xmapped  <- as.integer(parsedBamFile$xmapped) }
  if (!is.integer(parsedBamFile$xmismatch)) { parsedBamFile$xmismatch  <- as.integer(parsedBamFile$xmismatch) }
  if (!is.integer(parsedBamFile$xins)) { parsedBamFile$xins  <- as.integer(parsedBamFile$xins) }
  if (!is.integer(parsedBamFile$xdel)) { parsedBamFile$xdel  <- as.integer(parsedBamFile$xdel) }
  
  if (!is.numeric(parsedBamFile$mapq)) { parsedBamFile$mapq  <- as.numeric(parsedBamFile$mapq) }
  if (!is.numeric(parsedBamFile$readQ)) { parsedBamFile$readQ  <- as.numeric(parsedBamFile$readQ) }
  if (!is.numeric(parsedBamFile$readLen)) { parsedBamFile$readLen  <- as.numeric(parsedBamFile$readLen) }
  if (!is.numeric(parsedBamFile$meanCov)) { parsedBamFile$meanCov  <- as.numeric(parsedBamFile$meanCov) }
  if (!is.numeric(parsedBamFile$xcov)) { parsedBamFile$xcov <- as.numeric(parsedBamFile$xcov) }
  if (!is.numeric(parsedBamFile$xmapq)) { parsedBamFile$xmapq <- as.numeric(parsedBamFile$xmapq) }
  if (!is.numeric(parsedBamFile$xwidth)) { parsedBamFile$xwidth <- as.numeric(parsedBamFile$xwidth) }

  return(parsedBamFile)
}


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


# source(file.path("Scripts", "BamParser.R"))
# parsedBamFile <- parseBamFile(window.size=25000, force=TRUE)