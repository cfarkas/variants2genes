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
3) Call germline and somatic variants in both samples using Strelka2 variant caller.
4) Intersect (using --invert flag) strelka germline variants with bcftools filtered variants. 

- The output from these steps will output case-linked variants and correspondent genes. 
- Case-linked somatic variants and case-linked INDELs were be also reported.
```

# Installation:

### Option 1: Installing dependences via anaconda (recommended)
- requires miniconda, python2.7 and/or python>=3. To install miniconda, see: https://docs.conda.io/en/latest/miniconda.html
```
git clone https://github.com/cfarkas/variants2genes.git        # clone repository
cd variants2genes                                              # enter repository
conda config --add channels bioconda                           # add bioconda channel (if you haven't already done so)
conda env update --file environment.yml                        # install required programs
conda activate variants2genes                                  # load environment
bash makefile                                                  # make  & install
```
- Optionally (requires sudo privileges)
```
sudo cp ./bin/* /usr/local/bin/
```
After these steps, a conda enviroment called annotate_my_genomes can be managed as follows:
```
# To activate this environment, use
#
#     $ conda activate variants2genes
#
# To deactivate an active environment, use
#
#     $ conda deactivate
```

### Option 2: Without using conda, program by program:

- see detailed installation steps in our wiki here: https://github.com/cfarkas/annotate_my_genomes.wiki.git

# Usage:

After installation, provide:
- sorted Control and Case bam files (i.e.: WT and KO, respectively)
- genome assembly with correspondent GTF file
- number of processors.

Then, if variants2genes binary is copied in /usr/local/bin/, execute as follows:
```
variants2genes /path/to/WT.sorted.bam /path/to/KO.sorted.bam /path/to/mm10.fa /path/to/mm10.gtf 20
```
Otherwise, provide full path to this binary, located in ```variants2genes/bin/``` 

## Example: Collect haplotypes from RNA-seq data:
- As an example, we will analyze haplotypes from an RNA-seq data taken from SALL2 wild type and knockout mice, presenting germline variants linked to Chromosome 14, see: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5504-9. With the pipeline, we will obtain these linked variants to knockout mice, not present in the wild-type counterpart. The correspondent illumina reads will be downloaded and aligned against mm10 genome (mus musculus version 10). As outputs, the pipeline will take BAM file names until a point is encountered (i.e. for SRR8267474.sorted.bam ==> SRR8267474) so distinctive BAM file names are desired in the pipeline. After installation, inside variants2genes folder execute the following steps:

```
# Inside variants2genes folder
mkdir SALL2_WT_vs_KO

cd SALL2_WT_vs_KO

#######################
### Pipeline Starts ###
#######################

## STEP 1: Download reference genome from UCSC and correspondent GTF file. Then, build HISAT2 index. 
../bin/genome-download mm10
hisat2-build mm10.fa mm10_hisat2

## STEP 2: Download and align SALL2 Wild type and Knockout reads using SRA accessions, using 25 threads.

# WT sample
prefetch -O ./ SRR8267474 
fastq-dump --gzip SRR8267474
hisat2 -x mm10_hisat2 -p 25 -U SRR8267474.fastq.gz | samtools view -bSh > WT.bam

# KO sample
prefetch -O ./ SRR8267458 
fastq-dump --gzip SRR8267458
hisat2 -x mm10_hisat2 -p 25 -U SRR8267458.fastq.gz | samtools view -bSh > KO.bam

## STEP 3: sort bam samples using 25 threads
samtools sort -o WT.sorted.bam WT.bam -@ 25
samtools sort -o KO.sorted.bam KO.bam -@ 25

## STEP 4 (optional, but recommended): Use plot-variants to inspect genome-wide variants in every sample (check graph.pdf)
../bin/plot-variants WT.sorted.bam KO.sorted.bam mm10.fa ../R_scripts/bam_coverage_mouse.R

## STEP 5: Run variants2genes.sh script to collect KO-linked variants and correspondent genes with variants (using 20 threads)
../bin/variants2genes WT.sorted.bam KO.sorted.bam mm10.fa mm10.gtf 20
```
Check KO sub-folder with output files. From this example, two chr12 and 766 chr14 KO-linked germline variants were discovered.  

## Employing user-provided genome and/or GTF files:

Important: If users have their own genome and/or annotation file, their can use it in the pipeline, if desired. Their must edit STEP 1 and STEP 5. We will run the example using "my_genome.fa" and "my_annotation.gtf" instead of mm10.fa and mm10.gtf as follows:

```
## STEP 1: Build HISAT2 index of my_genome.fa. 
hisat2-build my_genome.fa my_genome_hisat2

## Same STEP 2-4

## STEP 5: Run variants2genes.sh script to collect KO-linked variants and correspondent genes with variants (using 20 threads)
./variants2genes WT.sorted.bam KO.sorted.bam my_genome.fa my_annotation.gtf 20
```

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc. 
