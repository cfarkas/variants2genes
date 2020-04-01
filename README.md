# variants2genes
Obtaining case-linked variants and correspondent genes (from control/case experiments) in BASH/R enviroment

![github_variants2genes](https://user-images.githubusercontent.com/7016350/77459123-d7d06d80-6dc4-11ea-8d21-54a6e7ca9c4b.png)


## Pipeline Outline

variants2genes are a set of bash scripts that address (based on several well-known genomic tools) case associated variants (with correspondent genes) 
between two matched samples from RNA-seq/WES/WGS data. The pipeline generates genome-wide plots of coverage between the two samples as initial inspection, and detect alleles present in the case sample (but not substantially present in the control) by calling variants and intersecting the correspondent genes using intersecting the output from bcftools and strelka somatic variant calling: https://github.com/Illumina/strelka. This resource could be useful in the Normal/Tumor comparison analysis, haplotype analysis, characterization of substrains, among other scenarios.
The pipeline is implemented in the BASH/R enviroment and is available for several organism models such as Human, Mouse, Rat and Chicken.

Pipeline Outline:

```
1) Call variants in Control and Case samples using bcftools mpileup.
2) Filter these variants using vcflib and bedtools (mainly to correct coverage artifacts)
3) Call germline and somatic variants in both samples using Strelka2 small variant caller.
4) Intersect (using --invert flag) strelka germline variants with bcftools filtered variants. 

- The output from these steps will output case-linked variants and correspondent genes. 
- Case-linked somatic variants and case-linked INDELs were be also reported.
```

## Installation requirements:
### Obtaining and installing R (>=3.2.0)
See https://cloud.r-project.org/ for R installation in linux/ubuntu. R version 3.2.3 comes from default in Ubuntu 16.04 LTS but users with older Ubuntu distributions must upgrade R. A way accomplish this can be the following:
```
# Removing R from system
sudo apt-get remove r-base-core

# Editing sources.list 
sudo su
echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Installing R (version 3.6.0)
sudo apt update; sudo apt install r-base
exit
R
```

The following R packages are required for the script usages:

dplyr >=0.7.4       
gridExtra >=2.3     
reshape2 >=1.4.2    
ggplot2 >=2.2.1 

These packages can be installed in R (macOS/ubuntu) by opening R and typying:
>install.packages("dplyr")<br/>install.packages("gridExtra")<br/>install.packages("reshape2")<br/>install.packages("ggplot2")
<br/>

To check installation of these packages, open R and type:

>library(dplyr)<br/>library(gridExtra)<br/>library(reshape2)<br/>library(ggplot2)<br/>

### Obtaining and installing vcflib:
Clone vcflib folder in current directory:
```
git clone --recursive git://github.com/vcflib/vcflib.git

# If this line not work, try:
git config --global url.https://github.com/.insteadOf git://github.com/
git clone --recursive git://github.com/vcflib/vcflib.git

#Enter vcflib directory and make
cd vcflib
make   # Needs CMake compiler, with sudo privileges do: sudo apt-get install cmake

#After make, binaries and scripts can be copied in /usr/local/bin with sudo. In vcflib/ directory:

sudo cp scripts/* /usr/local/bin/
sudo cp bin/* /usr/local/bin/

#To check vcflib scripts, type vcf in terminal followed by TAB and display all posibilities
```

### Obtaining and Installing BEDTools
Complete instructions can be found in https://bedtools.readthedocs.io/en/latest/content/installation.html. Users with privileges can accomplish with sudo: 
```
sudo apt-get install bedtools
```

### Obtaining and installing up-to-date SAMtools, bcftools and htslib (version 1.9)
Old samtools version will not work. Users needs to install version up to date of these three packages. Users can first install htslib v1.9 and then samtools with bcftools v1.9, respectively. For downloading these packages, see http://www.htslib.org/download/). The latter can be accomplish by downloading the three packages, decompressing it, and doing the following:
```
cd htslib-1.9    # and similarly for bcftools and samtools
sudo ./configure --prefix=/usr/local/bin
sudo make
sudo make install
# this step is only for samtools and bcftools...
sudo cp samtools /usr/local/bin/
```
Then in a terminal type
>samtools<br>bcftools

to check 1.9 versions (using htslib v1.9)

### Obtaining tabix (version >= 1.7, to decompress .vfc.gz files from bcftools output).
```
sudo apt install tabix
```

### Obtaining and installing BamTools
Complete instructions can be found in https://github.com/pezmaster31/bamtools/wiki/Building-and-installing. Users with privileges can accomplish with sudo: 

>sudo apt install bamtools

### Obtaining and installing Subread (for using featurecounts)
Complete instructions can be found in http://subread.sourceforge.net/. Users with privileges can accomplish with sudo: 

>sudo apt-get install subread

### Obtaining SRA toolkit from ncbi (for downloading reads from GEO datasets).
```
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz
gunzip sratoolkit.2.9.6-ubuntu64.tar.gz
tar -xvf sratoolkit.2.9.6-ubuntu64.tar
sudo cp sratoolkit.2.9.6-ubuntu64/bin/fastq-dump /usr/local/bin/

```
### Obtaining HISAT2 aligner (for aligning RNA-seq data).
```
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip
unzip hisat2-2.0.4-Linux_x86_64.zip
sudo cp hisat2-2.0.4/hisat2* /usr/local/bin/
```

# Usage:
## Collect haplotypes from RNA-seq data:
- As an example, we will analyze haplotypes from a pooled RNA-seq data taken from brain sections at 4 and 7 days of development. The correspondent illumina reads were aligned against galGal6 genome (gallus gallus version 6). The RNA sequencing comes from an unespecified gallus gallus substrain(s) without reference genome, and we need to dissect haplotypes (most likely germline variants) in both samples for variant characterization, when comparing within days. This will dissect alleles that varies in expression due genetic variation rather develpment. We will employ as target.bam the 7-day RNA sequencing, and the 4-day RNA sequencing as reference.bam. With these two bam files, the whole pipeline can be runned as follows:

```
git clone https://github.com/cfarkas/variants2genes
cd variants2genes
mkdir galGal6_analysis

# Copying all bash scripts to galGal6_analysis
cp bash_scripts/* ./galGal6_analysis/

# Copying relevant R script to galGal6_analysis
cp ./R_scripts/bam_coverage_chicken.R ./galGal6_analysis/

# Place reference and target RNA-seq bam file (illumina technology) into galGal6_analysis folder. IMPORTANT: bam files have to be named with a single word after .bam prefix. In this case we will name 4-day RNA-seq as "control.bam" and 7-day RNA-seq as "case.bam" 

cp /some_directory/control.bam ./galGal6_analysis/
cp /some_directory/case.bam ./galGal6_analysis/

#######################
### Pipeline Starts ###
#######################

## STEP 1: Inside galGal6_analysis folder, download reference genome from UCSC and correspondent GTF file.
cd galGal6_analysis/
bash genome_download.sh galGal6

## STEP 2: Use sort_bam.sh script to sort bam samples using 40 threads
bash sort_bam.sh control.bam case.bam 40

## STEP 3: Use plot-coverage.sh script to inspect genome-wide coverage (for "good" sequencing, check graph.pdf)
bash plot-coverage.sh control.sorted.bam case.sorted.bam bam_coverage_chicken.R 

## STEP 4: Run variants2genes.sh script to collect Case-linked variants and correspondent genes with variants (using 40 threads)
bash variants2genes.sh control.sorted.bam case.sorted.bam galGal6.fa galGal6.gtf 40

# All done. Check target sub-folder in ./galGal6_analysis with output files.
```

## Employing user-provided genome and/or GTF files:

Important: If users have their own genome and/or annotation file, their can use it in the pipeline, if desired. Their must skip STEP 1 and continue next steps. As example, we will use "my_genome.fa" and "final_annotated.gtf" instead of galGal6.fa and galGal6.gtf in STEP 1:

```
## STEP 2: Use sort_bam.sh script to sort bam samples using 40 threads
bash sort_bam.sh control.bam case.bam 40

## STEP 3: Use plot-coverage.sh script to inspect genome-wide coverage (for "good" sequencing, check graph.pdf)
bash plot-coverage.sh control.sorted.bam case.sorted.bam bam_coverage_chicken.R 

## STEP 4: Run variants2genes.sh script to collect Case-linked variants and correspondent genes with variants (using 40 threads)
bash variants2genes.sh control.sorted.bam case.sorted.bam my_genome.fa final_annotated.gtf 40

```
