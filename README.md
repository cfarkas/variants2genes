# variants2genes
Obtaining case-linked variants and correspondent genes (from control/case experiments)

![github_variants2genes](https://user-images.githubusercontent.com/7016350/77459123-d7d06d80-6dc4-11ea-8d21-54a6e7ca9c4b.png)

## Pipeline Outline

variants2genes is a pipeline that address (based on several well-known genomic tools) case associated variants (with correspondent genes) 
between two matched samples from RNA-seq/WES/WGS data. The pipeline generates genome-wide plots of coverage between the two samples as initial inspection, and detect alleles present in the case sample (but not substantially present in the control) by calling variants and intersecting the correspondent genes using intersecting the output from bcftools and strelka somatic variant calling: https://github.com/Illumina/strelka. This resource could be useful in the Normal/Tumor comparison analysis, haplotype analysis, characterization of substrains, among other scenarios.
The pipeline requieres BASH and R enviroment and is available for several organism models such as Human, Mouse, Rat and Chicken.

Pipeline Outline:

```
1) Picard and GATK4 best practices preprocessing of BAM files : https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
2) Call variants in Control and Case samples using bcftools mpileup, freebayes and gatk HaplotypeCaller.
3) Filter these variants using vcflib and bedtools (mainly to correct coverage artifacts)
4) Call germline and somatic variants in both samples using Strelka2 variant caller.
5) Intersect (using --invert flag) strelka somatic variants with bcftools/gatk HaplotypeCaller filtered variants.
6) join all variants for varlociraptor
7) Obtain somatic variants associated with case sample (FDR<=0.01) via varlociraptor filtering
8) Obtain germline variants shared between control and case samples (FDR<=0.01) via varlociraptor filtering
9) Obtain case-associated germline variants and corresponding genes harboring these variants (VCF and GTF format, respectively) by merging bcftools and gatk HaplotypeCaller filtered variants
```

# Installation:

### Option 1: Installing dependences via anaconda (recommended). Tested in Ubuntu 16.04, 20.04 and 22.04
- requires miniconda, python2.7 and/or python>=3. To install miniconda, see: https://docs.conda.io/en/latest/miniconda.html
```
git clone https://github.com/cfarkas/variants2genes.git        # clone repository
cd variants2genes                                              # enter repository
conda config --add channels bioconda                           # add bioconda channel (if you haven't already done so)
conda create --name variants2genes picard=2.18.7               # create environment
conda activate variants2genes                                  # activate environment

# Install packages
conda install -c conda-forge -y parallel 
conda install -c bioconda -y vcflib
conda install -c bioconda -y bedtools
conda install -c bioconda -y tabix
conda install -c bioconda -y sra-tools
conda install -c bioconda -y ncbi-ngs-sdk
conda install -c bioconda -y hisat2
conda install -c bioconda -y minimap2
conda install -c bioconda -y fastp
conda install -c conda-forge -y coreutils
conda install -c anaconda -y gawk
conda install -c conda-forge -y sed
conda install -c bioconda -y gatk4=4.0.5.1
conda install -c bioconda -y freebayes                                         
conda install -c anaconda -y pandas
conda install -c bioconda -y varlociraptor 
pip install cyvcf2
pip install notebook

bash makefile                                                  # make  & install
```
- Optionally (requires sudo privileges)
```
sudo cp ./bin/* /usr/local/bin/
```

Also install (not through conda):

- ```SAMtools``` and ```bcftools```, including  ```htslib```. To install it, see here: https://github.com/cfarkas/variants2genes/wiki#obtaining-and-installing-up-to-date-samtools-bcftools-and-htslib-latest-version115-march-24-2022
- ```Python2```. Some Ubuntu distros comes Python2 by default. In ubuntu 20.04 can be installed as follows: ```sudo apt install python2```

Optionally, install: 
- R packages ```ggplot2```, ```reshape2```, ```dplyr``` and ```gridextra```. See here: https://github.com/cfarkas/variants2genes/wiki#obtaining-and-installing-r-320


After these steps, a conda enviroment called variants2genes can be managed as follows:
```
# To activate this environment, use
#
#     $ conda activate variants2genes
#
# To deactivate an active environment, use
#
#     $ conda deactivate
```

#### Notes: 

- By activating variants2genes enviroment, all binaries in the variants2genes repository can be executed.

- Uninstall environment as follows: 
```
conda remove --name variants2genes --all
```

### Option 2: Without using conda, program by program:

- see detailed installation steps in our wiki here: https://github.com/cfarkas/annotate_my_genomes.wiki.git

# Usage:

Check binaries in ```variants2genes/bin/```. ```variants2genes``` pipeline will requiere:
```
    -a  File or path to Control bam file (sorted and indexed)
    -b  File or path to Case bam file (sorted and indexed)
    -g  Reference genome (in fasta format)
    -r  Gene annotation (in GTF format)
    -s  known SNPs for base recalibration (in VCF format)
    -t  Number of threads for processing (integer)
```
As example, for mouse RNA-seq data (mm10 genome), execute as follows:
```
# load environment
conda activate variants2genes

# Download mm10.fa and mm10.gtf files in place
genome-download mm10   

# Execute the pipeline using 20 threads for processing
variants2genes -a /path/to/WT.sorted.bam -b /path/to/KO.sorted.bam -g /path/to/mm10.fa -r /path/to/mm10.gtf -s dbSNP.raw.vcf -t 20
```
- After these steps, a folder named ```variants2genes_$DATE_OF_EXECUTION``` where ```$DATE_OF_EXECUTION = Day:Month:Year_Hour:Minute:Sec```, will contain the results. 
- To obtain mm10 known snps for mouse, please visit here: https://usegalaxy.org/u/carlosfarkas/h/dbsnpvcffiles
- To obtain mm39 known snps for mouse, please visit here:: https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/
- See an example in our wiki page here: https://github.com/cfarkas/variants2genes/wiki#example-collect-ko-linked-variants-from-rna-seq-data

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc. 
