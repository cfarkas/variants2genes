#!/bin/bash
{

control=${1}
case=${2}
ref=${3}
GTF=${4}
threads=${5}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {GTF} {Threads}"
  echo ""
  echo "This script will call variants using SAMtools/bcftools in Control and Case bam files to obtain case-ligated variants."
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Reference: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "GTF: PATH where the gene annotation (in GTF format) is located. If the GTF file is located in the working folder, just specify the name."
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {GTF} {Threads}"
  echo ""
  echo "This script will call variants using SAMtools/bcftools in Control and Case bam files to obtain case-ligated variants."
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Reference: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "GTF: PATH where the gene annotation (in GTF format) is located. If the GTF file is located in the working folder, just specify the name."
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {GTF} {Threads}"
  echo ""
  echo "This script will call variants using SAMtools/bcftools in Control and Case bam files to obtain case-ligated variants."
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Reference: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "GTF: PATH where the gene annotation (in GTF format) is located. If the GTF file is located in the working folder, just specify the name."
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {GTF} {Threads}"
  echo ""
  echo "This script will call variants using SAMtools/bcftools in Control and Case bam files to obtain case-ligated variants."
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Reference: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "GTF: PATH where the gene annotation (in GTF format) is located. If the GTF file is located in the working folder, just specify the name."
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {Control Bam file} {Case Bam file} {Reference} {GTF} {Threads}"; exit 1; }

if [ $# -ne 5 ]; then
  echo 1>&2 "Usage: bash ./`basename $0` {Control Bam file} {Case Bam file} {Reference} {GTF} {Threads}"
  exit 3
fi
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

### File name definitions
control_name=$(echo "${1}" | awk -F'[.]' '{print $1}')
case_name=$(echo "${2}" | awk -F'[.]' '{print $1}')

### Variant Calling
echo "................................................................................."
echo ""
echo "Performing Variant Calling with SAMtools and bcftools (see: http://samtools.github.io/bcftools/):"
echo ""
echo "The output directory will be the following:"
echo ${dir1}
begin=`date +%s`
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${ref} --threads ${threads} -Ou ${1}| bcftools call -mv -Ov -o ${control_name}.vcf
echo "done with Control Bam file. Continue with Case bam file..."
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${ref} --threads ${threads} -Ou ${2}| bcftools call -mv -Ov -o ${case_name}.vcf
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
echo "Filtering with bcftools control and case vcf files..."
bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 1' ${control_name}.vcf > Control_initial_filter.vcf
bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3' ${case_name}.vcf > ${case_name}.bcftools.vcf
echo "Done"
echo ""
echo "Selecting variants in case VCF not present in control VCF:"
vcfintersect -i Control_initial_filter.vcf ${case_name}.bcftools.vcf -r ${ref} --invert > case_variants.vcf
echo "Done"
echo ""
echo "Case VCF file: Filtering with vcflib ---> QUAL > 30"
vcffilter -f "QUAL > 30" case_variants.vcf > case_variants.QUAL1.filter.vcf
echo "First filter done"
echo ""
echo "Case VCF file: Filtering with vcflib ---> DP > 8"
vcffilter -f "DP > 8" case_variants.QUAL1.filter.vcf > case_variants.QUAL2.filter.vcf
echo "Second filter done"
echo ""
echo "Intersecting case variants in the ranges of control bam file:"
bamToBed -i ${1} > Control.bed
mergeBed -i Control.bed > Control.merged.bed
multiBamCov -bams ${1} -bed Control.merged.bed > control_counts
awk '{ if ($4 > 7) { print } }' control_counts > filter_merged.bed
vcfintersect -b filter_merged.bed case_variants.QUAL2.filter.vcf > Case.filtered.vcf
echo "Done"
echo "Initially filtered Case-associated variants are named Case.filtered.vcf"
rm Control.bed filter_merged.bed Control.merged.bed control_counts
rm case_variants*
rm *.filter1.vcf
rm Control_initial_filter.vcf
rm *bcftools.vcf
echo ""

echo "Performing Somatic Variant Calling with strelka v2.9.2:"
echo ""
echo "for documentation, please see: https://github.com/Illumina/strelka"
echo ""
echo "downloading strelka binary from github repository"
wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
tar xvjf strelka-2.9.2.centos6_x86_64.tar.bz2
echo ""
echo "run demo to check successful installation"
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaSomaticWorkflowDemo.bash
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaGermlineWorkflowDemo.bash
echo "demo run was succesfull"
echo "Running on control and case samples: Collecting Germline variants:"
echo ""
# configuration
begin=`date +%s`
./strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
    --bam ${1} \
    --bam ${2} \
    --referenceFasta ${3} \
    --runDir strelka_germline
# execution on a single local machine with n parallel jobs
strelka_germline/runWorkflow.py -m local -j ${5}
echo "Variant Calling done"
end=`date +%s`
elapsed=`expr $end - $begin`
echo "................................................................................."
cp ./strelka_germline/results/variants/variants.vcf.gz ./strelka_germline_variants.vcf.gz
bgzip -d strelka_germline_variants.vcf.gz
grep "#" strelka_germline_variants.vcf > strelka_germline_variants_header.vcf
grep "PASS" strelka_germline_variants.vcf > strelka_germline_variants_PASS.vcf
cat strelka_germline_variants_header.vcf strelka_germline_variants_PASS.vcf > strelka_germline_variants.filtered.vcf
rm strelka_germline_variants_header.vcf strelka_germline_variants_PASS.vcf
echo "Fitered variants are called strelka_germline_variants.filtered.vcf"
echo ""
echo "Continue with Somatic Variant Calling"
# configuration
begin=`date +%s`
./strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${1} \
    --tumorBam ${2} \
    --referenceFasta ${3} \
    --runDir strelka_somatic
# execution on a single local machine with n parallel jobs
strelka_somatic/runWorkflow.py -m local -j ${5}
echo "Variant Calling done"
end=`date +%s`
elapsed=`expr $end - $begin`
echo "................................................................................."
cp ./strelka_somatic/results/variants/somatic.snvs.vcf.gz ./strelka_somatic_variants.vcf.gz
bgzip -d strelka_somatic_variants.vcf.gz
grep "#" strelka_somatic_variants.vcf > strelka_somatic_variants_header.vcf
grep "PASS" strelka_somatic_variants.vcf > strelka_somatic_variants_PASS.vcf
cat strelka_somatic_variants_header.vcf strelka_somatic_variants_PASS.vcf > strelka_somatic_variants.filtered.vcf
rm strelka_somatic_variants_header.vcf strelka_somatic_variants_PASS.vcf
echo ""
echo "Filtered Germline and Somatic variants are located in working directory"
echo ""
echo "Filtering Case.filtered.vcf variants file with strelka outputs..."
vcfintersect -i strelka_germline_variants.filtered.vcf Case.filtered.vcf -r ${ref} > Case.filtered.strelka.vcf
echo "Done"

### Annotating variants and obtaining gene list
echo ""
echo "Annotating variants and the associated gene list"
bedtools intersect -a ${GTF} -b Case.filtered.strelka.vcf > Case.filtered.strelka.gtf
perl -lne 'print "@m" if @m=(/((?:gene_id)\s+\S+)/g);' Case.filtered.strelka.gtf > genes.with.variants.tabular
awk '!a[$0]++' genes.with.variants.tabular > gene_id_identifiers.tab
rm genes.with.variants.tabular
sed -i 's/gene_id //g' gene_id_identifiers.tab
sed -i 's/";//g' gene_id_identifiers.tab
sed -i 's/"//g' gene_id_identifiers.tab  
mv gene_id_identifiers.tab genes_with_variants.tabular
echo "done"

### Obtaining gene counts in both samples
echo ""
echo "Obtaining gene counts in both samples"
featureCounts -a ${GTF} -o gene_counts.tabular -T ${threads} ${1} ${2}
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
mv Case.filtered.strelka.gtf genes_with_variants.tabular Case.filtered.strelka.vcf count_table.tabular ${case_name}_quantification.tabular *.filtered.vcf ./${case_name}
echo "The following files are located in the the ./${case_name} folder"
echo ""
echo "(1) Case.filtered.strelka.gtf" 
echo "(2) genes_with_variants.tabular"    
echo "(3) Case.filtered.strelka.vcf"              
echo "(4) count_table.tabular"
echo "(5) ${case_name}_quantification.tabular"
echo "(6) strelka_germline_variants.filtered.vcf"
echo "(7) strelka_somatic_variants.filtered.vcf"
echo "(8) strelka_somatic_indels.filtered.vcf"
echo ""
echo "Corresponding to:"
echo ""
echo "(1): Case associated variants in GTF format"
echo "(2): List of genes with variants in tabular format"
echo "(3): Case-associated variants in VCF format"
echo "(4): Count table (gene_id) of Control and Case bam files"
echo "(5): List of genes with variants plus gene quantification in Control and Case bam files, respectively"
echo "(6): Strelka germline variants"
echo "(7): Strelka somatic variants associated with case bam file"
echo "(7): Strelka somatic indels associated with case bam file"
echo "....................................................................................................."

#
} | tee logfile
#
