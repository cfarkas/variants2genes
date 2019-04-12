# variants2genes
Obtaining case-associated variants and correspondent genes (from control/case experiments) in BASH/R enviroment

## Pipeline Outline

variants2genes are a set of bash scripts that address (based on several well-known genomic tools) case associated variants (with correspondent genes) 
between two matched samples from WES/WGS data (and in some cases, RNA-seq data). The pipeline generates genome-wide plots of coverage between the two samples as initial inspection, and detect alleles present in the case sample (but not in the control) by calling variants and intersecting the correspondent genes. This resource could be useful in the Normal/Tumor comparison analysis, among other scenarios.
The pipeline is implemented in the BASH/R enviroment and is available for several organism models such as Human, Mouse and Rat.

## Preeliminars:
### Obtaining and installing R
See https://www.r-project.org/ for R installation. R comes by default in Ubuntu but could be outdated. Users can also accomplish via conda:

>conda install -c r r 

The following R packages are required for the script usages (see https://www.r-project.org/ for R installation)

dplyr >=0.7.4       
gridExtra >=2.3     
reshape2 >=1.4.2    
ggplot2 >=2.2.1 

These packages can be installed in R (macOS/ubuntu) by opening R and typying:
>install.packages("dplyr")<br/>install.packages("gridExtra")<br/>install.packages("reshape2")<br/>install.packages("ggplot2")
<br/>

To check installation of these packages, open R and type:

>library(dplyr)<br/>library(gridExtra)<br/>library(reshape2)<br/>library(ggplot2)<br/>

### Obtaining and installing Freebayes:
Clone Freebayes folder in current directory: 
>git clone --recursive git://github.com/ekg/freebayes.git

Enter Freebayes directory and make:
>cd freebayes/<br/> sudo make install

Copy scripts to /usr/local/bin:
>sudo cp ./scripts/* /usr/local/bin/

To check installation, type in terminal:
>freebayes<br/>bamleftalign

### Obtaining and installing vcflib:
Clone vcflib folder in current directory:
>git clone --recursive git://github.com/vcflib/vcflib.git

Enter vcflib directory and make
>cd vcflib/<br/>make   

After make, binaries and scripts can be copied in /usr/local/bin with sudo. In vcflib/ directory:

>sudo cp scripts/* /usr/local/bin/<br/>sudo cp bin/* /usr/local/bin/

To check vcflib scripts, type vcf in terminal followed by TAB and display all posibilities

### Obtaining and Installing BEDTools
Complete instructions can be found in https://bedtools.readthedocs.io/en/latest/content/installation.html. Users with privileges can accomplish with sudo: 

>sudo apt-get install bedtools

### Obtaining and installing SAMtools
Complete instructions can be found in http://www.htslib.org/. Users with privileges can accomplish with sudo: 

>sudo apt-get install samtools

### Obtaining and installing BamTools
Complete instructions can be found in https://github.com/pezmaster31/bamtools/wiki/Building-and-installing. Users with privileges can accomplish with sudo: 

>sudo apt install bamtools 

### Obtaining and installing bcftools
Complete instructions can be found in https://samtools.github.io/bcftools/. Users with privileges can accomplish with sudo: 

>sudo apt install bcftools

### Obtaining and installing Subread
Complete instructions can be found in http://subread.sourceforge.net/. Users with privileges can accomplish with sudo: 

>sudo apt-get install subread
