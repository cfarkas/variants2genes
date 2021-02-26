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


# Installation:

Just three steps in a terminal:

```
git clone https://github.com/cfarkas/variants2genes.git
cd variants2genes
bash makefile
```

# Usage:
## Collect haplotypes from RNA-seq data:
- As an example, we will analyze haplotypes from an RNA-seq data taken from SALL2 wild type and knockout mice, presenting germline variants linked to Chromosome 14, see: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5504-9. With the pipeline, we will obtain these linked variants to knockout mice, not present in the wild-type counterpart. The correspondent illumina reads will be downloaded and aligned against mm10 genome (mus musculus version 10). It is important to denominate every bam file with a name without points in the middle (e.g.: WT.bam vs KO.bam or Control.bam vs Case.bam and NOT Control.1.bam) since the pipeline will fail with complex names and points in it. After installation, inside variants2genes folder execute the following steps:

```
# Inside variants2genes folder
mkdir SALL2_WT_vs_KO

# Copying all binaries to SALL2_WT_vs_KO
cp ./bin/* ./SALL2_WT_vs_KO/

# Copying relevant R script to SALL2_WT_vs_KO
cp ./R_scripts/bam_coverage_mouse.R ./SALL2_WT_vs_KO/

cd SALL2_WT_vs_KO

#######################
### Pipeline Starts ###
#######################

## STEP1: Download reference genome from UCSC and correspondent GTF file. Then, build HISAT2 index. 
./genomeDownload mm10
hisat2-build mm10.fa mm10_hisat2

## STEP2: Align SALL2 Wild type and Knockout reads using SRA accessions. Give to bam files simple names (WT.bam and KO.bam) 

hisat2 -x mm10_hisat2 -p 25 --sra-acc SRR8267474,SRR8267475,SRR8267476,SRR8267477 | samtools view -bSh > WT.bam
hisat2 -x mm10_hisat2 -p 25 --sra-acc SRR8267458,SRR8267459,SRR8267460,SRR8267461 | samtools view -bSh > KO.bam

## STEP 3: Use sort_bam.sh script to sort bam samples using 40 threads
./sortBam WT.bam KO.bam 25

## STEP 4 (optional, but recommended): Use plotVariants to inspect genome-wide variants in every sample (check graph.pdf)
./plotVariants WT.sorted.bam KO.sorted.bam mm10.fa bam_coverage_mouse.R 

## STEP 5: Run variants2genes.sh script to collect KO-linked variants and correspondent genes with variants (using 20 threads)
./variants2genes WT.sorted.bam KO.sorted.bam mm10.fa mm10.gtf 20

# All done. Check KO sub-folder with output files.
```
From this example, two chr12 and 766 chr14 KO-linked germline variants were discovered.  

## Employing user-provided genome and/or GTF files:

Important: If users have their own genome and/or annotation file, their can use it in the pipeline, if desired. Their must edit STEP1 and STEP5. We will run the example using "my_genome.fa" and "final_annotated.gtf" instead of mm10.fa and mm10.gtf as follows:

```
## STEP1: Download reference genome from UCSC and correspondent GTF file. Then, build HISAT2 index. 
hisat2-build my_genome.fa my_genome_hisat2

## STEP2: Align SALL2 Wild type and Knockout reads using SRA accessions. Give to bam files simple names (WT.bam and KO.bam) 

hisat2 -x my_genome_hisat2 -p 25 --sra-acc SRR8267474,SRR8267475,SRR8267476,SRR8267477 | samtools view -bSh > WT.bam
hisat2 -x my_genome_hisat2 -p 25 --sra-acc SRR8267458,SRR8267459,SRR8267460,SRR8267461 | samtools view -bSh > KO.bam

## STEP 3: Use sort_bam.sh script to sort bam samples using 40 threads
./sortBam WT.bam KO.bam 25

## STEP 4: Use plotVariants to inspect genome-wide variants in every sample (check graph.pdf)
./plotVariants WT.sorted.bam KO.sorted.bam mm10.fa bam_coverage_mouse.R 

## STEP 5: Run variants2genes.sh script to collect KO-linked variants and correspondent genes with variants (using 20 threads)
./variants2genes WT.sorted.bam KO.sorted.bam my_genome.fa final_annotated.gtf 20
```

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc. 
