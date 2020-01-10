# variants2genes
Obtaining case-linked variants and correspondent genes (from control/case experiments) in BASH/R enviroment

## Pipeline Outline

variants2genes are a set of bash scripts that address (based on several well-known genomic tools) case associated variants (with correspondent genes) 
between two matched samples from WES/WGS data (RNA-seq data can be also used). The pipeline generates genome-wide plots of coverage between the two samples as initial inspection, and detect alleles present in the case sample (but not substantially present in the control) by calling variants and intersecting the correspondent genes. This resource could be useful in the Normal/Tumor comparison analysis, haplotype analysis, characterization of substrains, among other scenarios.
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

# Quickstart:
## Collect haplotypes from RNA-seq data:
- In this example, we will analyze haplotypes from an RNA-seq data, aligned againts galGal6 genome (gallus gallus version 6). This RNA sequencing comes from an unespecified gallus gallus substrain (named in this example as aligned.sorted.bam), without reference genome. To call haplotypes in this sample, we will employ as reference a WGS illumina bam used in the assembly of galGal6 genome (sra accession SRR3954707) . The whole pipeline can be runned as follows:
```
git clone https://github.com/cfarkas/variants2genes
cd variants2genes
mkdir galGal6_analysis

# Copying all bash scripts to galGal6_analysis
cp bash_scripts/* ./galGal6_analysis/

# Copying relevant R script to galGal6_analysis
cp ./R_scripts/bam_coverage_chicken.R ./galGal6_analysis/

# Place reference (illumina WGS bam file) and target RNA-seq bam file into galGal6_analysis folder. IMPORTANT: bam files have to be named with a single word after .bam prefix. In this case we will name reference WGS sequencing as "reference.bam" and target RNA-seq sequencing as "target.bam" 
cp /some_directory/reference.bam ./galGal6_analysis/
cp /some_directory/target.bam ./galGal6_analysis/

# Download Reference genome from UCSC into galGal6_analysis folder
bash genome_download.sh galGal6

#######################
### Pipeline Starts ###
#######################

## STEP 1: Use sort_bam.sh script to sort bam samples using 40 threads
bash sort_bam.sh reference.bam target.bam 40

## STEP 2: Use plot-coverage.sh script to inspect genome-wide coverage (check graph.pdf)
bash plot-coverage.sh reference.sorted.bam target.sorted.bam bam_coverage_chicken.R 

## STEP 3: Run variants2genes.sh script to collect Case-linked variants and correspondent genes with variants (using 40 threads)
bash variants2genes.sh reference.sorted.bam target.sorted.bam galGal6.fa 40

# All done. Check target sub-folder in ./galGal6_analysis with output files.
```
