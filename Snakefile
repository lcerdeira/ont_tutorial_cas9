import os
import re
import pandas as pd
import sys
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

usableThreads = 8

# characterise and/or prepare the specified reference genome
reference_genome = config["reference_genome"]
# ReferenceGenome is a derivative of the defined location
ReferenceGenome = os.path.join("ReferenceData", os.path.basename(reference_genome))
ReferenceGenomeMMI = re.sub("(.zip)|(.bz2)|(.gz)", "", ReferenceGenome)+".mmi"
gzippedRef = False
UnzippedReferenceGenome = None
# accommodate some handling of compressed genome references
if (re.search("\.gz$", ReferenceGenome)):
  gzippedRef = True
  UnzippedReferenceGenome = re.sub("(\.zip$)|(\.bz2$)|(\.gz$)", "", ReferenceGenome)

# some fastq handling - a fastq may exist as a provided file that may be compressed or uncompressed; or may be a pointer
# to a folder that contains one or more fastq files that may themselves be compressed or uncompressed
# since these files could be quite massive I would prefer to avoid duplication, uncompress or other IO unncessary tasks ...
fastq = config["fastq"]
fastqTarget = config["study_name"]

# load targetIds from supplied bed file
df = pd.read_csv(config["target_regions"], sep="\t", header=None)
TARGETS = list(df[df.columns[3]])

# decide how we will be accessing the FASTQ files ...
fastqAccessParam = ""
if (os.path.isfile(fastq) and re.search("(.gz$)|(.gzip$)", fastq) is not None): # compressed input -file-
  #print("InputFastq is compressed file")
  fastqAccessParam = "gunzip -c "+fastq
elif (os.path.isfile(fastq) and re.search("(.gz$)|(.gzip$)", fastq) is None): # uncompressed fastq -file-
  #print("InputFastq is uncompressed file")
  fastqAccessParam = "cat "+fastq
elif (os.path.isdir(fastq)):
  files = os.listdir(fastq)
  fastqRx = re.compile("(\.fastq$)|(\.fq$)")
  fastqZRx = re.compile("(\.fastq.gz$)|(\.fq.gz$)")
  fastqNative = len(list(filter(fastqRx.search, files)))
  fastqCompre = len(list(filter(fastqZRx.search, files)))
  if (fastqNative>=fastqCompre):
    #print("InputFastq is set of uncompressed file*s*")
    fastqAccessParam = "cat "+ " ".join([fastq + os.path.sep + x for x in list(filter(fastqRx.search, files))])
  else:
    #print("InputFastq is set of compressed file*s*")
    fastqAccessParam = "gunzip -c "+ " ".join([fastq + os.path.sep + x for x in list(filter(fastqRx.search, files))])
#print(fastqAccessParam)


"""
This snakemake rule will download a reference genome from an open genome database etc
- assuming that this will be a fasta file that may or may not be compressed ...
"""
rule download_reference_genome:
  input:
    HTTP.remote(reference_genome, keep_local=True)
  output:
    ReferenceGenome
  run:
    shell("mv {input} {output}")


"""
A simple method to quickly unpack a gzipped reference; lazy
- this only makes sense for the IGV step in the workflow
"""
rule unzip_reference_genome:
  input:
    ReferenceGenome
  output:
    UnzippedReferenceGenome
  run:
    shell("gunzip -c {input} > {output}")



"""
create Minimap2 index - this has a minimal additional workflow cost; reap benefits if the workflow is run multiple times
- this does increase diskspace requirements - but ... you shouldn't map without a load of space?
"""
# build minimap2 index
rule Minimap2Index: 
  input:
    genome = UnzippedReferenceGenome if gzippedRef else ReferenceGenome
  output:
    index = ReferenceGenomeMMI
  threads: usableThreads
  shell:
    """
    minimap2 -t {threads} -d {output.index} {input.genome}
    """


"""
run Minimap2 on the provided sequence file(s)
"""
rule Minimap2:
  input:
    fastq = fastq,
    ref = ReferenceGenomeMMI
  output:
    bam = "Analysis/Minimap2/"+fastqTarget+".bam",
    ubam = "Analysis/Minimap2/"+fastqTarget+".unmapped.bam"
  params:
    rg = config["pipeline"],
    fap = fastqAccessParam
  threads: usableThreads
  shell:
    """{params.fap} | minimap2 -2 -a -x map-ont --MD -R '@RG\\tID:{params.rg}\\tSM:{params.rg}' -t {threads} {input.ref} - | samtools view -@ {threads} -F 0x4 -O BAM -U {output.ubam} | samtools sort -@ {threads} > {output.bam}"""


# build BAM index
rule SortedBamIndex: 
  input:
    bam = "Analysis/Minimap2/"+fastqTarget+".bam"
  output:
    index = "Analysis/Minimap2/"+fastqTarget+".bam.bai"
  threads: usableThreads
  shell:
    "samtools index -@ {threads} -b {input.bam}"


# render quality values from unmapped reads 
rule RenderUnmappedQvals:
  input:
    ubam = "Analysis/Minimap2/"+fastqTarget+".unmapped.bam"
  output:
    uqual = "Analysis/Minimap2/"+fastqTarget+".unmapped.quals"
  threads: 8
  shell:
    "samtools view -@ {threads} -O sam {input.ubam} | awk '{{print $11}}' > {output.uqual}"


# perform the R preprocess 
rule Rpreprocess:
  input:
    bam = "Analysis/Minimap2/"+fastqTarget+".bam",
    bai = "Analysis/Minimap2/"+fastqTarget+".bam.bai",
    uqual = "Analysis/Minimap2/"+fastqTarget+".unmapped.quals"
  output:
    rout = "Analysis/R/"+fastqTarget+"_mapping_results.Rdata",
    targets = expand("Analysis/OnTarget/"+fastqTarget+".{target}.mappedreads", target=TARGETS),
    offtargets = "Analysis/OffTarget/"+fastqTarget+".OffTarget.mappedreads"
  shell:
    "Rscript Static/Source/harvest.R"
    
# filter the starting fastq for ontarget reads
rule onTargetReadDump:
  input:
    rin = "Analysis/R/"+fastqTarget+"_mapping_results.Rdata", # not actually used; for method chaining only
    readIds = "Analysis/OnTarget/"+fastqTarget+".{target}.mappedreads"
  output:
    fastq = "Analysis/OnTarget/"+fastqTarget+".{target}.fastq"
  params:
    fap = fastqAccessParam
  shell:
    "{params.fap} | seqtk subseq - {input.readIds} > {output.fastq}"

# filter the starting fastq for the offtarget reads
rule offTargetReadDump:
  input:
    rin = "Analysis/R/"+fastqTarget+"_mapping_results.Rdata", # not actually used; for method chaining only
    readIds = "Analysis/OffTarget/"+fastqTarget+".OffTarget.mappedreads"
  output:
    fastq = "Analysis/OffTarget/"+fastqTarget+".OffTarget.fastq"
  params:
    fap = fastqAccessParam
  shell:
    "{params.fap} | seqtk subseq - {input.readIds} > {output.fastq}"


rule all:
  input:
    expand("Analysis/OnTarget/"+fastqTarget+".{target}.fastq", target=TARGETS),
    "Analysis/OffTarget/"+fastqTarget+".OffTarget.fastq"
