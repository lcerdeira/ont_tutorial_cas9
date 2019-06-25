# Statement of tutorial objectives

The aim of this tutorial is to demonstrate a workflow for mapping long DNA sequence reads to a reference genome, and to evaluate the performance of a Cas9 based target enrichment strategy. This workflow is suitable for Oxford Nanopore fastq sequence collections and requires a reference genome and a BED file of target coordinates.

The tutorial is packaged with example data, and the workflow can be reproduced to address questions such as

* How many sequence reads map to the reference genome?
* What is the depth of coverage for reads that map to pre-defined target regions?
* What is the background depth of coverage for the non-targetted genomic regions?
* Is there evidence for off-target enrichment of genomic regions?
* How has the Cas9 based target-enrichment worked for **`HTT`**, our gene of interest?

Editing the workflow's configuration file, **`config.yaml`**, will allow the analyses to be run using different DNA sequence collections, reference genomes, and with different BED files that define the regions of interest.

## What will this tutorial produce?

* A rich HTML format report containing summary statistics and figures highlighting performance of the enrichment protocol
* **`Microsoft Excel`** format files containing coordinates and summary statistics for on-target regions
* **`Fastq`** format file containing reads that correspond to each of the target regions specified
* **`Microsoft Excel`** format files containing summary statistics for sequence regions defined as showing off-target enrichment
* **`Fastq`** format sequence file containing reads that were classified as off-target enrichment
* **`Coordinates`** and instructions for reviewing candidate genomic regions with **`IGV`** the Integrated Genomics Viewer


## Methods utilised include: 

* **`conda`** for management of bioinformatics software installations
* **`snakemake`** for managing the bioinformatics workflow
* **`minimap2`** for mapping sequence reads to reference genome
* **`samtools`** for SAM/BAM handling and mapping statistics
* **`RSamtools`** and **`GenomicAlignments`**; R software for parsing BAM files
* **`seqtk`** for writing out the subseq of sequence reads that map to the target region
* **`IGV`** for visualising mapping characteristics at specific genomic regions

## The computational requirements include: 

* Computer running Linux (Centos7, Ubuntu 18_10, Fedora 29, macOS)
* Multiple CPU cores are ideal; a 4 core CPU at minimum would be recommended 
* At least 16 Gb RAM - this is sufficient for mapping against the human genome and can report a full MinION flowcell worth of sequence data. The packaged dataset uses just human chromosome 4 and 8Gb RAM is sufficient
* At least 15 Gb spare disk space for analysis and indices
* Runtime with provided example data - approximately 45 minutes

\pagebreak

# Software installation

1. Most software dependencies are managed though **`conda`**. Install as described at  <br> [https://conda.io/docs/install/quick.html](https://conda.io/docs/install/quick.html). You will need to accept the license agreement during installation and we recommend that you allow the Conda installer to prepend its path to your `.bashrc` file when asked.
```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    bash
```
2. Download Nanopore tutorial & example files into a folder named `ont_tutorial_cas9`. This tutorial requires the **`git-lfs`** large file support capabilities, which should be installed through **`Conda`** first
```
    conda install -c conda-forge git-lfs
    git lfs install
    git clone https://github.com/nanoporetech/ont_tutorial_cas9.git ont_tutorial_cas9
```
3. Change your working directory into the new `ont_tutorial_cas9` folder 
```
    cd ont_tutorial_cas9
```
4. Install Conda software dependencies
```
    conda env create --name ont_tutorial_cas9 --file environment.yaml
```
5. Initialise the Conda environment 
```
    source activate ont_tutorial_cas9
    R CMD javareconf > /dev/null 2>&1
```


\pagebreak

# Introduction

Targeted sequencing strategies provide a cost-effective way of sequencing regions of interest to high coverage. Unlike other enrichment methods nanopore sequencing does not require amplification of any sort which allows us to target:

* Long gene targets that are not amenable to long-range PCR (>30 kb) either in a single pass or using up to 100 target sites in a single assay (known as a tiling approach)
* Regions with methylation patterns or modifications (that can be recovered using software such as TOMBO)
* Regions that are highly repetitive

[Oxford Nanopore Technologies](https://nanoporetech.com) provides a [Cas-mediated PCR-free enrichment protocol](https://community.nanoporetech.com/protocols/Cas-mediated-PCR-free-enrich/) which uses an enrichment strategy based on the design of CRISPR RNA (crRNA) probe sequences that may flank or tile-across one or more target regions. The crRNAs program the Cas9 protein to bind and cleave DNA at sites that match the crRNA sequence. This Cas9-mediated cut of the DNA, and the production of a newly-exposed and “deprotected” DNA end is the basis of the enrichment protocol.

For an equivalent amount of sequenced library a Cas-mediated DNA enrichment will provide a higher coverage for the targeted regions than sequencing the whole genome alone.

The *Cas-mediated PCR-free enrichment protocol* recommends a strategy for the design of CRISPR RNA (crRNA) probes and provides recommendations for the quality of the purified DNA. These recommendations aim to maximise the on-target sequence recovery by reducing both the DNA background and amounts of off-target DNA. When performing optimally the protocol will delivers up to a couple of Gb of data from a MinION flow cell from 48 hours of sequencing and provide 100s-1000x coverage across target regions. This tutorial has been prepared to help assess performance of a Cas-mediated enrichment study and to help identify target regions and DNA preparation steps that may require further optimisation.


There are five goals for this tutorial:

* To introduce a literate framework for analysing Oxford Nanopore *Cas-mediated PCR-free enriched DNA sequence* data
* To utilise best data-management practices
* To map sequence reads to the reference genome
* To identify and report the sequence reads that map to the defined target regions and off-target regions
* To describe the relative depletion of the background genome


# Getting started and best practices

This tutorial requires a computer workstation running a Linux operating system. The workflow described has been tested using **`Fedora 29`**, **`Centos 7`** and **`Ubuntu 18_04`**. This tutorial has been prepared in the **`Rmarkdown`** file format. This utilises *markdown* (an easy-to-write plain text format as used in many Wiki systems) - see @R-rmarkdown for more information about **`rmarkdown`**. The document template contains chunks of embedded **`R code`** that are dynamically executed during the report preparation. 

The described analytical workflow makes extensive use of the **`conda`** package management and the **`snakemake`** workflow software. These software packages and the functionality of **`Rmarkdown`** provide the source for a rich, reproducible and extensible tutorial document.

The workflow contained within this Tutorial performs an authentic bioinformatics analysis and using the human chromosome 4 as a reference sequence. There are some considerations in terms of memory and processor requirement. Indexing the whole human genome for sequence read mapping using **`minimap2`** will use at least **`18 Gb`** of memory. The minimal recommended hardware setup for this tutorial is a 4 threaded computer with at least 8 Gb of RAM and 10 Gb of storage space. 

There are few dependencies that need to be installed at the system level prior to running the tutorial. The **`conda`** package management software will coordinate the installation of the required bioinformatics software and their dependencies in user space - this is dependent on a robust internet connection.

As a best practice this tutorial will separate primary DNA sequence data (the base-called fastq files) from the **`Rmarkdown`** source and the genome reference data. The analysis results and figures will again be placed in a separate working directory. The required layout for the primary data is shown in the figure below. This minimal structure will be prepared over the next few sections of this tutorial. The DNA sequences and mapping BED file must be placed within a folder called **`RawData`** and the reference genome and annotation files must be placed in a folder named **`ReferenceData`**. All results will be placed in a folder named **`Analysis`** and different sub-folders will be created for different steps in the workflow. The **`Static`** folder contains accessory methods, texts, bibliography and graphics.


![](Static/Images/FolderLayout.png) 

# Experimental setup

The first required step for performing a meta-analysis of a *Cas-mediated PCR-free enrichment protocol* based sequencing study is to define the experimental design. 

The example data included with this tutorial describes a Cas-mediated enrichment experiment that targets the HTT gene. The enriched sequence library has been sequenced using a single MinION flowcell on a GridION sequencing device. To manage the sizes of the datasets downloaded, the example dataset provided has been filtered to select for only the sequence reads on **`Chromosome 4`**. The example data is therefore a **synthetic dataset** - but no other enrichment or modification of the sequences has been performed.  

The design for the tutorial is defined within a YAML format configuration file (**`config.yaml`**). The tutorial's file is shown in the figure below.


![](Static/Images/ExperimentalDesign.png) 




* **`pipeline`** identifies the workflow that this configuration file belongs to
* **`study_name`** is a label used to identify the analysis and the label is used to name the files produced during the analysis
* **`reference_genome`** refers to the genome against which the sequence reads will be mapped. In this tutorial a URL is provided and the **`Snakemake`** workflow will download the corresponding file.
* **`target_regions`** points to a BED format file that describes the genomic coordinates for each of the targets being assessed. The BED file format is a tab-delimited file with un-named columns describing, in order, the chromosome, start position, end position and target name.
* **`fastq`** is a pointer to either a single fastq file (may be gzipped) or a folder of fastq files. These are the sequences that will be mapped to the **`reference_genome`** and assessed for mapping to the targets defined in **`target_regions`**
* **`gstride`**, **`target_proximity`** and **`offtarget_level`** are parameters used to control the analysis. For offtarget and background assessment, the genome is split into windows **`gstride`** nucleotides in length and mean coverage is assessed. Regions with a mean-coverage of > **`offtarget_level X`** the mean background level are defined as being off-target enrichment regions. **`target_proximity`** defines the window up- and down-stream of the **`target_regions`** that are excluded from the background calculations. This reduces overall mean background by excluding ontarget reads that extend beyond the target-region.


\newpage

# Snakemake

This tutorial for the assessment of target enrichment from DNA sequence data makes use of **`snakemake`** (@snakemake2012). Snakemake is a workflow management system implemented in Python. The aim of the Snakemake tool is to enable reproducible and scalable data analyses. The workflow produced within this document should be portable between laptop resources, computer servers and other larger scale IT deployments. The Snakemake workflow additionally defines the sets of required software (and software versions where appropriate) and will automate the installation and deployment of this software through the **conda** package management system.

The **`snakemake`** workflow will call methods that include **`minimap2`** @minimap22018 and **`samtools`** @samtools2009. The planned workflow is shown in the figure below. 


![](Static/Images/dag.png) 

The precise commands within the **`Snakefile`** include

* download the specified reference genome
* use **`minimap2`** to index the reference genome
* map DNA sequence reads against the reference genome index using **`minimap2`**
* convert **`minimap2`** output (**`SAM`**) into a sorted **`BAM`** format using **`samtools`** and filter out the unmapped reads
* prepare summary mapping statistics using **`Rsamtools`** and **`GenomicAlignments`** from the provided **`R`** script
* filter for the on-target sequence reads using **`seqtk`** (@seqtkurl)

# Run the snakemake workflow file


\fontsize{8}{12}
```
# just type snakemake to run the workflow
# Methods run inside the workflow can utilise multiple compute cores for parallel performance
# Use the -j flag to specify the number of cores to use; the example below requests 4 cores

snakemake -j 4 all
```
\fontsize{10}{14}


\pagebreak
 

# Prepare the analysis report

The **`Rmarkdown`** script can be run using the **`knit`** dialog in the **`Rstudio`** software. 

The document can also be rendered from the command line with the following command. This command is  run automatically during the Snakemake workflow.

\fontsize{8}{12}
```
R --slave -e 'rmarkdown::render("ont_tutorial_cas9.Rmd", "html_document")'
```
\fontsize{10}{14}

