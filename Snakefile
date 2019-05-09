import os
import re
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

# characterise and/or prepare the specified reference genome
reference_genome = config["reference_genome"]
# ReferenceGenome is a derivative of the defined location
ReferenceGenomeIncoming = os.path.join("ReferenceData", os.path.basename(reference_genome))
ReferenceGenome = re.sub("(.zip)|(.bz2)|(.gz)", "", os.path.join("Analysis", "ReferenceGenome", os.path.basename(reference_genome)))
ReferenceGenomeMMI = ReferenceGenome+".mmi"

# some fastq handling - a fastq may exist as a provided file that may be compressed or uncompressed; or may be a pointer
# to a folder that contains one or more fastq files that may themselves be compressed or uncompressed
# since these files could be quite massive I would prefer to avoid duplication, uncompress or other IO unncessary tasks ...
fastq = config["fastq"]
fastqTarget = config["study_name"]

"""
This snakemake rule will download a reference genome from an open genome database etc
- assuming that this will be a fasta file that may or may not be compressed ...
"""
rule download_reference_genome:
  input:
    HTTP.remote(reference_genome, keep_local=True)
  output:
    ReferenceGenomeIncoming
  run:
    shell("mv {input} {output}")


"""
The supplied genome reference may be gzip / bzip2 / ... compressed - snakemake has problems
with uncompressing files into the same folder - therefore we will either copy or uncompress
the file into a location within the Analysis folder - this is not required for Minimap2 or
NGMLR
"""
rule unpack_reference_genome:
  input:
    ReferenceGenomeIncoming
  output:
    ReferenceGenome
  run:
    if (re.search("\.gz$", ReferenceGenomeIncoming)):
      shell("gunzip -c {input} > {output}")
    elif (re.search("\.bz2$", ReferenceGenomeIncoming)):
      shell("bunzip2 --keep -d {input}")
    else:
      shell("cp {input} {output}")



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
    rg = config["pipeline"]
  run:
    source = None
    # this code doesn't look very easy to read - there is considerable duplication
    # I am not sure how we could for example, force and evaluation of a variable (might contain --cat--) within a shell command
    if (os.path.isfile(fastq) and re.search("(.gz$)|(.gzip$)", fastq) is not None): # compressed input -file-
      shell("""gunzip -c {fastq} | minimap2 -2 -a -x map-ont --MD -R '@RG\tID:{params.rg}\tSM:{params.rg}' -t 8 {input.ref} - | samtools view -@ 4 -F 0x4 -O BAM -U {output.ubam} | samtools sort -@ 4 > {output.bam}""")
    elif (os.path.isfile(fastq) and re.search("(.gz$)|(.gzip$)", fastq) is None): # uncompressed fastq -file-
      shell("""cat {fastq} | minimap2 -2 -a -x map-ont --MD -R '@RG\tID:{params.rg}\tSM:{params.rg}' -t 8 {input.ref} - | samtools view -@ 4 -F 0x4 -O BAM -U {output.ubam} | samtools sort -@ 4 > {output.bam}""")
    elif (os.path.isdir(fastq)):
      files = os.listdir(fastq)
      fastqRx = re.compile("(\.fastq$)|(\.fq$)")
      fastqZRx = re.compile("(\.fastq.gz$)|(\.fq.gz$)")
      fastqNative = len(list(filter(fastqRx.search, files)))
      fastqCompre = len(list(filter(fastqZRx.search, files)))
      if (fastqNative>=fastqCompre):
        shell("cat "+ " ".join([fastq + os.path.sep + x for x in list(filter(fastqRx.search, files))])+""" | minimap2 -2 -a -x map-ont --MD -R '@RG\\tID:{params.rg}\\tSM:{params.rg}' -t 8 {input.ref} - | samtools view -@ 4 -F 0x4 -O BAM -U {output.ubam} | samtools sort -@ 4 > {output.bam}""")
      else:
        shell("gunzip -c "+ " ".join([fastq + os.path.sep + x for x in list(filter(fastqZRx.search, files))])+""" | minimap2 -2 -a -x map-ont --MD -R '@RG\tID:{params.rg}\tSM:{params.rg}' -t 8 {input.ref} - | samtools view -@ 4 -F 0x4 -O BAM -U {output.ubam} | samtools sort -@ 4 > {output.bam}""")      


# build minimap2 index
rule SortedBamIndex: 
  input:
    bam = "Analysis/Minimap2/"+fastqTarget+".bam"
  output:
    index = "Analysis/Minimap2/"+fastqTarget+".bam.bai"
  shell:
    "samtools index -b {input.bam}"



rule all:
  input:
    ReferenceGenome,
    "Analysis/Minimap2/"+fastqTarget+".bam.bai"
