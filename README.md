# variants2genes
Obtaining case-linked variants and correspondent genes (from control/case experiments) in BASH/R enviroment

## Pipeline Outline

variants2genes are a set of bash scripts that address (based on several well-known genomic tools) case associated variants (with correspondent genes) 
between two matched samples from WES/WGS data (and in some cases, RNA-seq data). The pipeline generates genome-wide plots of coverage between the two samples as initial inspection, and detect alleles present in the case sample (but not substantially present in the control) by calling variants and intersecting the correspondent genes. This resource could be useful in the Normal/Tumor comparison analysis, among other scenarios.
The pipeline is implemented in the BASH/R enviroment and is available for several organism models such as Human, Mouse and Rat.

## Preeliminars:
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
make   

#After make, binaries and scripts can be copied in /usr/local/bin with sudo. In vcflib/ directory:

sudo cp scripts/* /usr/local/bin/
sudo cp bin/* /usr/local/bin/

#To check vcflib scripts, type vcf in terminal followed by TAB and display all posibilities
```

### Obtaining and Installing BEDTools
Complete instructions can be found in https://bedtools.readthedocs.io/en/latest/content/installation.html. Users with privileges can accomplish with sudo: 

>sudo apt-get install bedtools

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

### Obtaining and installing BamTools
Complete instructions can be found in https://github.com/pezmaster31/bamtools/wiki/Building-and-installing. Users with privileges can accomplish with sudo: 

>sudo apt install bamtools

### Obtaining and installing Subread (for using featurecounts)
Complete instructions can be found in http://subread.sourceforge.net/. Users with privileges can accomplish with sudo: 

>sudo apt-get install subread

### Obtaining SRA toolkit from ncbi (for downloading reads for Quick-start).
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

# Quick start (mouse RNA-seq):
```
git clone https://github.com/cfarkas/variants2genes
cd variants2genes

####################
### Preeliminars ###
####################

# Downloading test fastq files in test folder.
mkdir test
cd test/
fastq-dump -Z SRR8267474 > WT.fastq
fastq-dump -Z SRR8267458 > KO1.fastq
cd ..

# Coping bash scripts and bam_coverage_mouse.R script to test folder. 
cp bash_scripts/* ./test/
cp ./R_scripts/bam_coverage_mouse.R ./test/

# Download mm10 genome, index it and download mm10 GTF annotation file (using genome_download.sh script).
cd test/
bash genome_download.sh mm10    # for mouse genome mm10 build.
mkdir hisat2_index_mm10

# Build genome index using 40 threads (-p parameter).
hisat2-build mm10.fa -p 40 ./hisat2_index_mm10/mm10_index

# Align reads to reference genome using 40 threads (-p parameter).
hisat2 -x ./hisat2_index_mm10/mm10_index -p 40 -U WT.fastq | samtools view -bS - > WT.bam
hisat2 -x ./hisat2_index_mm10/mm10_index -p 40 -U KO1.fastq | samtools view -bS - > KO1.bam

#######################
### Pipeline Starts ###
#######################

## STEP 1: Use sort_bam.sh script to sort bam samples using 40 threads
bash sort_bam.sh WT.bam KO1.bam 40

## STEP 2: Use plot-coverage.sh script to inspect genome-wide coverage (check graph.pdf)
bash plot-coverage.sh WT.sorted.bam KO1.sorted.bam bam_coverage_mouse.R 

## STEP 3: Run variants2genes.sh script to collect Case-linked variants and correspondent genes with variants (using 40 threads)
bash variants2genes.sh WT.sorted.bam KO1.sorted.bam mm10.fa 40

# All done. Check KO1 sub-folder in ./test with output files.
```
