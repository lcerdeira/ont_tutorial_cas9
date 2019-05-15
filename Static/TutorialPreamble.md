# Statement of tutorial objectives

The aim of this tutorial is to demonstrate a workflow for mapping long DNA sequence reads to a reference genome, and to evaluate the performance of a Cas9 based target enrichment strategy. This workflow is suitable for Oxford Nanopore fastq sequence collections but requires a reference genome and a BED file of target coordinates.

The tutorial is packaged with example data, so that the workflow can be replicated to address questions such as

* How many sequence reads map to the reference genome?
* What is the depth of coverage for reads that map to pre-defined target regions?
* What is the background depth of coverage for the non-targetted genomic regions?
* Is there evidence for off-target enrichment of genomic regions?
* How has the Cas9 based target-enrichment worked for **`HTT`**, our gene of interest?

Editing of the workflow's configuration file, **`config.yaml`**, will allow the workflow to be run with different DNA sequence collections, reference genomes, and with different BED files that define the regions of interest.

## Methods utilised include: 

* **`conda`** for management of bioinformatics software installations
* **`snakemake`** for managing the bioinformatics workflow
* **`minimap2`** for mapping sequence reads to reference genome
* **`samtools`** for SAM/BAM handling and mapping statistics
* **`RSamtools`** and **`GenomicAlignments`**; R software for parsing BAM files
* **`seqtk`** for writing out the subseq of sequence reads that map to the target region

## The computational requirements include: 

* Computer running Linux (Centos7, Ubuntu 18_10, Fedora 29)
* At least 24 Gb RAM - this is sufficient for a single thread and can accommodate a MinION flowcell worth of sequence data
* At least 15 Gb spare disk space for analysis and indices
* Runtime with provided example data - approximately 20 minutes

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
```


\pagebreak

# Introduction

PLEASE CAN I HAVE A COUPLE OF PARAGRAPHS THAT DESCRIBE

1. THE REASON WHY WE HAVE A CAS9 PROTOCOL
2. WHAT WE EXPECT A CAS9 ENRICHMENT TO DELIVER
3. SOME OVERLAP OF THE PROTOCOL WITH PRODUCTS AND STUFF ON THE STORE



There are five goals for this tutorial:

* To introduce a literate framework for analysing Oxford Nanopore DNA sequence data prepared using the MinION, GridION or PromethION
* To utilise best data-management practices
* To provide basic DNA sequence QC metrics, enabling review and consideration of the starting experimental data
* To map sequence reads to the reference genome and to identify the assess target regions, target-proximal regions and the background-genome for patterns of read mapping
* To assess target regions for depth-of-coverage and to explore strandedness of mapping


# Getting started and best practices

This tutorial requires a computer workstation running a Linux operating system. The workflow described has been tested using **`Fedora 29`**, **`Centos 7`** and **`Ubuntu 18_04`**. This tutorial has been prepared in the **`Rmarkdown`** file format. This utilises *markdown* (an easy-to-write plain text format as used in many Wiki systems) - see @R-rmarkdown for more information about **`rmarkdown`**. The document template contains chunks of embedded **`R code`** that are dynamically executed during the report preparation. 

The described analytical workflow makes extensive use of the **`conda`** package management and the **`snakemake`** workflow software. These software packages and the functionality of **`Rmarkdown`** provide the source for a rich, reproducible and extensible tutorial document.

The workflow contained within this Tutorial performs an authentic bioinformatics analysis and using the whole human genome as a reference sequence. There are some considerations in terms of memory and processor requirement. Indexing the whole human genome for sequence read mapping using **`minimap2`** for example will use at least **`18 Gb`** of memory. The minimal recommended hardware setup for this tutorial is therefore a 4 threaded computer with at least 24 Gb of RAM and 15 Gb of storage space. 

There are few dependencies that need to be installed at the system level prior to running the tutorial. The **`conda`** package management software will coordinate the installation of the required bioinformatics software and their dependencies in user space - this is dependent on a robust internet connection.

As a best practice this tutorial will separate primary DNA sequence data (the base-called fastq files) from the **`Rmarkdown`** source and the genome reference data. The analysis results and figures will again be placed in a separate working directory. The required layout for the primary data is shown in the figure below. This minimal structure will be prepared over the next few sections of this tutorial. The DNA sequences must be placed within a folder called **`RawData`** and the reference genome and annotation files must be placed in a folder named **`ReferenceData`**.


![](Static/Images/FolderLayout.png) 

# Experimental setup

The first required step for performing a  differential sequence analysis involves collation of information on the sequence collection, the reference genome and the target regions that should be enriched within the sequence collection.

![](Static/Images/ExperimentalDesign.png) 

The example data included with this tutorial describes a Cas9 experiment that has been used to enrich the HTT gene target. The enriched sequence library has been sequenced using a MinION flowcell on a GridION sequencing machine. To keep the sizes of the datasets downloaded more manageable, the example dataset provided has been filtered to select for only the sequence reads on **`Chromosome 4`**. Similarly, the **`snakemake`** workflow provided will download just the chromosome 4 reference sequence rather than the whole genome. The example data is therefore a **synthetic dataset** - but no other enrichment or modification of the sequences has been performed.  

This design is described in a configuration file named **`config.yaml`** - an example file has been provided with the tutorial. The content of this file is highlighted in the figure above. 

**`reference_genome`** refers to the genome against which the sequence reads will be mapped; **`genome_annotation`** refers to the gene annotations assigned to this genome sequence. In this tutorial a URL is provided for both and the **`Snakemake`** workflow will download the corresponding files. 

The only other parameters that should be considered include

* **`target_regions`**, 
* **`study`**,
* **`target_proximity`**,
* **`offtarget_level`**


\newpage

# Snakemake

This tutorial for the assessment of target enrichment from DNA sequence data makes use of **`snakemake`** (@snakemake2012). Snakemake is a workflow management system implemented in Python. The aim of the Snakemake tool is to enable reproducible and scalable data analyses. The workflow produced within this document should be portable between laptop resources, computer servers and other larger scale IT deployments. The Snakemake workflow additionally defines the sets of required software (and software versions where appropriate) and will automate the installation and deployment of this software through the **conda** package management system.

The **`snakemake`** workflow will call methods that include **`minimap2`** @minimap22018 and **`samtools`** @samtools2009. The planned workflow is shown in the figure below. The provided reference genome sequence will be indexed using **`minimap2`**, each sequence collection will be mapped to the index (again using **`minimap2`** with parameters tuned for the mapping of long reads whilst accommodating exon matches interspersed by introns) and summary statistics will be prepared using the **`samtools`** software. The remainder of the analysis will be performed in the **`R analysis`** described within the report.

![](Static/Images/dag.png) 

The precise commands within the **`Snakefile`** include

* download the specified reference genome
* download the specified genome annotations
* use **`minimap2`** to index the reference genome
* map DNA sequence reads against the reference genome index using **`minimap2`**
* convert **`minimap2`** output (**`SAM`**) into a sorted **`BAM`** format using **`samtools`**
* prepare summary mapping statistics using **`Rsamtools`** in a provided **`R`** script
* filter for the on-target sequence reads using **`seqtk`**

# Run the snakemake workflow file


\fontsize{8}{12}
```
# just type snakemake to run the workflow
# don't type <NPROC> but specify the number of processor cores available (e.g. 2 or 4)

snakemake -j <NPROC> all
```
\fontsize{10}{14}


\pagebreak
 

# Prepare the analysis report

The **`Rmarkdown`** script can be run usimg the **`knit`** dialog in the **`Rstudio`** software. 

The document can also be rendered from the command line with the following command. This command is also run by the Snakemake workflow

\fontsize{8}{12}
```
R --slave -e 'rmarkdown::render("ont_tutorial_cas9.Rmd", "html_document")'
```
\fontsize{10}{14}

