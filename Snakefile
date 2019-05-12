import os
import re
import pandas as pd
import sys
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

# characterise and/or prepare the specified reference genome
reference_genome = config["reference_genome"]
# ReferenceGenome is a derivative of the defined location
ReferenceGenome = os.path.join("ReferenceData", os.path.basename(reference_genome))
ReferenceGenomeMMI = re.sub("(.zip)|(.bz2)|(.gz)", "", ReferenceGenome)+".mmi"


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
  print("InputFastq is compressed file")
  fastqAccessParam = "gunzip -c "+fastq
elif (os.path.isfile(fastq) and re.search("(.gz$)|(.gzip$)", fastq) is None): # uncompressed fastq -file-
  print("InputFastq is uncompressed file")
  fastqAccessParam = "cat "+fastq
elif (os.path.isdir(fastq)):
  files = os.listdir(fastq)
  fastqRx = re.compile("(\.fastq$)|(\.fq$)")
  fastqZRx = re.compile("(\.fastq.gz$)|(\.fq.gz$)")
  fastqNative = len(list(filter(fastqRx.search, files)))
  fastqCompre = len(list(filter(fastqZRx.search, files)))
  if (fastqNative>=fastqCompre):
    print("InputFastq is set of uncompressed file*s*")
    fastqAccessParam = "cat "+ " ".join([fastq + os.path.sep + x for x in list(filter(fastqRx.search, files))])
  else:
    print("InputFastq is set of compressed file*s*")
    fastqAccessParam = "gunzip -c "+ " ".join([fastq + os.path.sep + x for x in list(filter(fastqRx.search, files))])
print(fastqAccessParam)


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
create Minimap2 index - this has a minimal additional workflow cost; reap benefits if the workflow is run multiple times
- this does increase diskspace requirements - but ... you shouldn't map without a load of space?
"""
# build minimap2 index
rule Minimap2Index: 
  input:
    genome = ReferenceGenome
  output:
    index = ReferenceGenomeMMI
  shell:
    "minimap2 -t 8 -d {output.index} {input.genome}"


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
  shell:
    """{params.fap} | minimap2 -2 -a -x map-ont --MD -R '@RG\\tID:{params.rg}\\tSM:{params.rg}' -t 8 {input.ref} - | samtools view -@ 4 -F 0x4 -O BAM -U {output.ubam} | samtools sort -@ 4 > {output.bam}"""


# build BAM index
rule SortedBamIndex: 
  input:
    bam = "Analysis/Minimap2/"+fastqTarget+".bam"
  output:
    index = "Analysis/Minimap2/"+fastqTarget+".bam.bai"
  shell:
    "samtools index -b {input.bam}"


# render quality values from unmapped reads 
rule RenderUnmappedQvals:
  input:
    ubam = "Analysis/Minimap2/"+fastqTarget+".unmapped.bam"
  output:
    uqual = "Analysis/Minimap2/"+fastqTarget+".unmapped.quals"
  shell:
    "samtools view -@ 5 -O sam {input.ubam} | awk '{{print $11}}' > {output.uqual}"


# perform the R preprocess 
rule Rpreprocess:
  input:
    bam = "Analysis/Minimap2/"+fastqTarget+".bam",
    bai = "Analysis/Minimap2/"+fastqTarget+".bam.bai",
    uqual = "Analysis/Minimap2/"+fastqTarget+".unmapped.quals"
  output:
    rout = "Analysis/R/"+fastqTarget+"_mapping_results.Rdata",
    targets = expand("Analysis/OnTarget/{target}.mappedreads", target=TARGETS),
    offtargets = "Analysis/OnTarget/OffTarget.mappedreads"
  shell:
    "Rscript harvest.R"
    
# filter the starting fastq for ontarget reads
rule onTargetReadDump:
  input:
    rin = "Analysis/R/"+fastqTarget+"_mapping_results.Rdata", # not actually used; for method chaining only
    readIds = "Analysis/OnTarget/{target}.mappedreads"
  output:
    fastq = "Analysis/OnTarget/{target}.fastq"
  params:
    fap = fastqAccessParam
  shell:
    "{params.fap} | seqtk subseq - {input.readIds} > {output.fastq}"

# filter the starting fastq for the offtarget reads
rule offTargetReadDump:
  input:
    rin = "Analysis/R/"+fastqTarget+"_mapping_results.Rdata", # not actually used; for method chaining only
    readIds = "Analysis/OnTarget/OffTarget.mappedreads"
  output:
    fastq = "Analysis/OnTarget/OffTarget.fastq"
  params:
    fap = fastqAccessParam
  shell:
    "{params.fap} | seqtk subseq - {input.readIds} > {output.fastq}"


# render the report ...
rule renderTheReport:
  input:
    rin = "Analysis/R/"+fastqTarget+"_mapping_results.Rdata",
    fastq = "Analysis/OnTarget/{target}.fastq",
    ofastq = "Analysis/OnTarget/OffTarget.fastq"
  output:
    "ont_tutorial_cas9.html"
  shell:
    "R --slave -e 'rmarkdown::render("ont_tutorial_cas9.Rmd", "html_document")'"


rule all:
  input:
    expand("Analysis/OnTarget/{target}.fastq", target=TARGETS),
    "Analysis/OnTarget/OffTarget.fastq",
    "ont_tutorial_cas9.html"
