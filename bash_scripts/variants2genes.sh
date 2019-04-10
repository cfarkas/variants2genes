#!/bin/bash

control=${1}
case=${2}
ref=${3}
threads=${4}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {Threads}"
  echo ""
  echo "This script will call variants using freebayes-parallel in Control and Case bam files to obtain case-ligated variants."
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Reference: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {Threads}"
  echo ""
  echo "This script will call variants using freebayes-parallel in Control and Case bam files to obtain case-ligated variants"
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Reference: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {Threads}"
  echo ""
  echo "This script will call variants using freebayes-parallel in Control and Case bam files to obtain case-ligated variants"
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Reference: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {Threads}"
  echo ""
  echo "This script will call variants using freebayes-parallel in Control and Case bam files to obtain case-ligated variants"
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Reference: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {Control Bam file} {Case Bam file} {Reference} {Threads}"; exit 1; }

if [ $# -ne 4 ]; then
  echo 1>&2 "Usage: bash ./`basename $0` {Control Bam file} {Case Bam file} {Reference} {Threads}"
  exit 3
fi
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

### File name definitions
control_name=$(echo "${1}" | awk -F'[.]' '{print $1}')
case_name=$(echo "${2}" | awk -F'[.]' '{print $1}')

### Variant Calling
echo "................................................................................."
echo ""
echo "Performing Variant Calling with Freebayes (see: https://github.com/ekg/freebayes):"
echo ""
echo "The output directory will be the following:"
echo ${dir1}
begin=`date +%s`
#freebayes-parallel <(fasta_generate_regions.py ${ref}.fai 100000) ${threads} -f ${ref} -b ${1} > ${control_name}.vcf
echo "done with Control Bam file. Continue with Case bam file..."
#freebayes-parallel <(fasta_generate_regions.py ${ref}.fai 100000) ${threads} -f ${ref} -b ${2} > ${case_name}.vcf
end=`date +%s`
elapsed=`expr $end - $begin`
echo ""
echo "Variant Calling done"
echo "................................................................................."
echo Time taken: $elapsed
echo ""
echo "Filtering and intersecting VCF files..."
echo ""

### Filtering and intersecting VCF files
echo "Initial filtering for control VCF file (QUAL > 1):"
vcffilter -f "QUAL > 1" ${control_name}.vcf > Control_initial_filter.vcf
echo "done"
echo "Selecting variants in case VCF not present in control VCF:"
vcfintersect -i Control_initial_filter.vcf ${case_name}.vcf -r ${ref} --invert > case_variants.vcf
echo "done"
echo "Filtering with QUAL > 30 case_variants.vcf:"
vcffilter -f "QUAL > 30" case_variants.vcf > case_variants.QUAL.vcf
echo "done"
echo "Filtering with DP > 5:"
vcffilter -f "DP > 5"  case_variants.QUAL.vcf > case_variants.QUAL.DP.vcf
echo "done"
echo "Intersecting case variants in the ranges of control bam file:"
bamToBed -i ${1} > Control.bed
mergeBed -i Control.bed > merged.bed
vcfintersect -b merged.bed case_variants.QUAL.DP.vcf > Case.filter.vcf
echo "done"
echo "Decomposing complex variants:"
vcfallelicprimitives -g Case.filter.vcf > Case.filtered.vcf
echo "done"
echo "Filtered Case-associated variants are named Case.filtered.vcf"
rm merged.bed
rm Control.bed
rm case_variants*
rm *.filter.vcf
rm Control_initial_filter.vcf

### Annotating variants and obtaining gene list
echo ""
echo "Annotating variants and the associated gene list"
genome_name=$(echo "${3}" | awk -F'[.]' '{print $1}')
bedtools intersect -a ${genome_name}.gtf -b Case.filtered.vcf > ${case_name}_annot_variants.gtf
cat ${case_name}_annot_variants.gtf | awk '{print $10}' > genes1.tabular
sed 's/"//' genes1.tabular > genes2.tabular
sed 's/";//' genes2.tabular > genes3.tabular
sort genes3.tabular | uniq -u > genes_with_variants.tabular
rm genes1.tabular genes2.tabular genes3.tabular
echo ""
echo "Done. ${case_name}_annot_variants.gtf contains the intersected variants with the reference genome"
echo "genes_with_variants.tabular is the list of genes with variants"

### Obtaining gene counts in both samples
echo ""
echo "Obtaining gene counts in both samples"
featureCounts -a ${genome_name}.gtf -o gene_counts.tabular -T ${threads} ${1} ${2}
cat gene_counts.tabular | awk '{print $1"\t"$7"\t"$8}' > count_table.tabular
rm gene_counts.tabular gene_counts.tabular.summary
echo "Done. count_table.tabular file contains the gene quantification across bam files"
echo ""
echo "Adding gene count detection to the list of genes with variants"
join <(sort genes_with_variants.tabular) <(sort count_table.tabular) > ${case_name}_quantification.tabular
echo ""
echo "All done. ${case_name}_quantification.tabular contains the list of genes with case associated-variants plus gene expression quantification in both samples"
echo ""
mkdir ${case_name} 
mv ${case_name}_annot_variants.gtf genes_with_variants.tabular ${control_name}.vcf ${case_name}.vcf Case.filtered.vcf count_table.tabular ${case_name}_quantification.tabular ./${case_name}
echo "The following files are located in the the ./${case_name} folder"
echo ""
echo "(1) ${case_name}_annot_variants.gtf" 
echo "(2) genes_with_variants.tabular"    
echo "(3) ${control_name}.vcf" 
echo "(4) ${case_name}.vcf"
echo "(5) Case.filtered.vcf"              
echo "(6) count_table.tabular"
echo "(7) ${case_name}_quantification.tabular"
echo ""
echo "Corresponding to:"
echo ""
echo "(1): Case associated variants in GTF format"
echo "(2): List of genes with variants in tabular format"
echo "(3): Raw Control VCF file"
echo "(4): Raw Case VCF file"
echo "(5): Case-associated variants in VCF format"
echo "(6): Count table (gene_id) of Control and Case bam files"
echo "(7): List of genes with variants plus gene quantification in Control and Case bam files, respectively"
echo "....................................................................................................."  





