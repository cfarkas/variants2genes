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

### Option 1: Installing dependences via anaconda (recommended)
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
pip install notebook                                           
pip install pandas
pip install cyvcf2                                       

bash makefile                                                  # make  & install
```
- Optionally (requires sudo privileges)
```
sudo cp ./bin/* /usr/local/bin/
```

Also install (not through conda):

- ```SAMtools``` and ```bcftools```. To install it, see here: https://github.com/cfarkas/variants2genes/wiki#obtaining-and-installing-up-to-date-samtools-bcftools-and-htslib-latest-version115-march-24-2022
- ```Varlociraptor``` (NOT through conda, please use ```cargo```). To install it, see here: https://varlociraptor.github.io/docs/installation/ and here: https://github.com/cfarkas/variants2genes/wiki#obtaining-rust-cargo-rust-package-and-install-velociraptor
- ```Python2```. Some Ubuntu distros comes Python2 by default. In ubuntu 20.04 can be installed as follows: ```sudo apt install python2```
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
variants2genes -a /path/to/WT.sorted.bam -b /path/to/KO.sorted.bam -g /path/to/mm10.fa -r /path/to/mm10.gtf -s mm10_dbSNP.raw.vcf -t 20
```
- After these steps, a folder named ```variants2genes_$DATE_OF_EXECUTION``` where ```$DATE_OF_EXECUTION = Day:Month:Year_Hour:Minute:Sec```, will contain the results. 
- To obtain known snps for each species, please visit here: https://usegalaxy.org/u/carlosfarkas/h/dbsnpvcffiles

## Example: Collect KO-linked variants from RNA-seq data:
- As an example, we will analyze haplotypes from an RNA-seq data taken from SALL2 wild type and knockout mice, presenting germline variants linked to Chromosome 14, see: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5504-9. With the pipeline, we will obtain these linked variants to knockout mice, not present in the wild-type counterpart. The correspondent illumina reads will be downloaded and aligned against mm10 genome (mus musculus version 10). 
- As outputs, the pipeline will take BAM file names until the ```.bam``` prefix is encountered
- From scratch, inside ```variants2genes``` folder:

#### Obtaining BAM files to be used as inputs:
```
# Inside variants2genes folder
mkdir SALL2_WT_vs_KO && cd SALL2_WT_vs_KO

## (1) Download reference genome from UCSC and correspondent GTF file. Then, build HISAT2 index. 
../bin/genome-download mm10
hisat2-build mm10.fa mm10_hisat2

## (2) Download and align SALL2 Wild type and Knockout reads with HISAT2, using 25 threads.
prefetch -O ./ SRR8267474 && fastq-dump --gzip SRR8267474            # WT sample
hisat2 -x mm10_hisat2 -p 25 -U SRR8267474.fastq.gz | samtools view -bSh > WT.bam
prefetch -O ./ SRR8267458 && fastq-dump --gzip SRR8267458            # KO sample
hisat2 -x mm10_hisat2 -p 25 -U SRR8267458.fastq.gz | samtools view -bSh > KO.bam

## (3) sort and index bam samples using 25 threads
samtools sort -o WT.sorted.bam WT.bam -@ 25 && samtools index WT.sorted.bam -@ 25
samtools sort -o KO.sorted.bam KO.bam -@ 25 && samtools index KO.sorted.bam -@ 25

## (4) (optional): Use plot-variants to inspect genome-wide variants in every sample (check graph.pdf)
../bin/plot-variants -a WT.sorted.bam -b KO.sorted.bam -g mm10.fa -p ../R_scripts/bam_coverage_mouse.R
```
#### Running the pipeline:
```
## (5) Run variants2genes.sh pipeline to collect KO-linked variants and correspondent genes with variants (using 20 threads)
wget -O mm10_dbSNP.raw.vcf https://usegalaxy.org/datasets/bbd44e69cb8906b509fb7398dabbcd16/display?to_ext=vcf   # download mm10 known-snps sites
../bin/variants2genes -a WT.sorted.bam -b KO.sorted.bam -g mm10.fa -r mm10.gtf -s mm10_dbSNP.raw.vcf -t 20
```
Inside ```variants2genes_$DATE_OF_EXECUTION```, check ```output_files``` sub-folder containing output files. From this example, two chr12 and 502 chr14 KO-linked germline variants were discovered.  

## Employing user-provided genome and/or GTF files:

Important: If users have their own genome and/or annotation file, their can use it in the pipeline, if desired. Their must edit STEP 1 and STEP 5. We will run the example using ```my_genome.fa``` and ```my_annotation.gtf``` instead of ```mm10.fa``` and ```mm10.gtf``` as follows:

#### Obtaining BAM files to be used as inputs:
```
## (1): Build HISAT2 index of my_genome.fa. 
hisat2-build my_genome.fa my_genome_hisat2

... Same Steps (2), (3) and (4) ...
```
#### Running the pipeline:
```
## (5) Run variants2genes.sh script to collect KO-linked variants and correspondent genes with variants (using 20 threads)
wget -O mm10_dbSNP.raw.vcf https://usegalaxy.org/datasets/bbd44e69cb8906b509fb7398dabbcd16/display?to_ext=vcf   # download mm10 known-snps sites
../bin/variants2genes -a WT.sorted.bam -b KO.sorted.bam -g my_genome.fa -r my_annotation.gtf -s mm10_dbSNP.raw.vcf -t 20
```

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc. 
