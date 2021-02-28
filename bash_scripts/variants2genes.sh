#!/bin/bash

set -e 

{

control=${1}
case=${2}
ref=${3}
GTF=${4}
threads=${5}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [Control Bam File] [Case Bam File] [Reference] [GTF] [Threads]"
  echo ""
  echo "This program will call variants using SAMtools/bcftools in Control and Case bam files to obtain case-ligated variants and will filter these variants by using BEDtools and Strelka somatic variant caller."
  echo ""
  echo "[Control Bam File]: File of path to Control bam file (sorted and indexed)"
  echo ""
  echo "[Case Bam File]: File of path to Case bam file (sorted and indexed)"
  echo ""
  echo "[Reference]: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[GTF]: PATH where the gene annotation (in GTF format) is located. If the GTF file is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {GTF} {Threads}"
  echo ""
  echo "This program will call variants using SAMtools/bcftools in Control and Case bam files to obtain case-ligated variants and will filter these variants by using BEDtools and Strelka somatic variant caller."
  echo ""
  echo "[Control Bam File]: File of path to Control bam file (sorted and indexed)"
  echo ""
  echo "[Case Bam File]: File of path to Case bam file (sorted and indexed)"
  echo ""
  echo "[Reference]: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[GTF]: PATH where the gene annotation (in GTF format) is located. If the GTF file is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {GTF} {Threads}"
  echo ""
  echo "This program will call variants using SAMtools/bcftools in Control and Case bam files to obtain case-ligated variants and will filter these variants by using BEDtools and Strelka somatic variant caller."
  echo ""
  echo "[Control Bam File]: File of path to Control bam file (sorted and indexed)"
  echo ""
  echo "[Case Bam File]: File of path to Case bam file (sorted and indexed)"
  echo ""
  echo "[Reference]: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[GTF]: PATH where the gene annotation (in GTF format) is located. If the GTF file is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` {Control Bam File} {Case Bam File} {Reference} {GTF} {Threads}"
  echo ""
  echo "This program will call variants using SAMtools/bcftools in Control and Case bam files to obtain case-ligated variants and will filter these variants by using BEDtools and Strelka somatic variant caller."
  echo ""
  echo "[Control Bam File]: File of path to Control bam file (sorted and indexed)"
  echo ""
  echo "[Case Bam File]: File of path to Case bam file (sorted and indexed)"
  echo ""
  echo "[Reference]: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[GTF]: PATH where the gene annotation (in GTF format) is located. If the GTF file is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [Control Bam file] [Case Bam file] [Reference] [GTF] [Threads]"; exit 1; }

if [ $# -ne 5 ]; then
  echo 1>&2 "Usage: ./`basename $0` [Control Bam file] [Case Bam file] [Reference] [GTF] [Threads]"
  exit 3
fi
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
begin=`date +%s`
#    .---------- constant part!
#    vvvv vvvv-- the code from above
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
### File name definitions
control_name=$(echo "${1}" | awk -F'[.]' '{print $1}')
case_name=$(echo "${2}" | awk -F'[.]' '{print $1}')
### Variant Calling
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Performing Variant Calling with bcftools (see: http://samtools.github.io/bcftools/):"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
echo ""
begin=`date +%s`
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${ref} --threads ${threads} -Ou ${1}| bcftools call -mv -Ov -o ${control_name}.vcf
echo "done with Control Bam file. Continue with Case bam file..."
echo ""
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${ref} --threads ${threads} -Ou ${2}| bcftools call -mv -Ov -o ${case_name}.vcf
end=`date +%s`
elapsed=`expr $end - $begin`
echo ""
printf "${CYAN}:::::::::::::::::::::\n"
echo "Variant Calling done"
echo ""
echo Time taken: $elapsed
printf "${CYAN}:::::::::::::::::::::${NC}\n"
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Filtering and intersecting VCF files"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::${NC}\n"
echo ""
### Filtering and intersecting VCF files
echo "Filtering with bcftools control and case vcf files..."
bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 1' ${control_name}.vcf > Control_initial_filter.vcf
bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3' ${case_name}.vcf > ${case_name}.bcftools.vcf
echo "Done"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Selecting variants in case VCF not present in control VCF:"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
vcfintersect -i Control_initial_filter.vcf ${case_name}.bcftools.vcf -r ${ref} --invert > case_variants.vcf
echo "Done"
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Case VCF file: Filtering with vcflib ---> QUAL > 30"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
vcffilter -f "QUAL > 30" case_variants.vcf > case_variants.QUAL1.filter.vcf
echo "First filter done"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Case VCF file: Filtering with vcflib ---> DP > 8"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
vcffilter -f "DP > 8" case_variants.QUAL1.filter.vcf > case_variants.QUAL2.filter.vcf
echo "Second filter done"
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Intersecting case variants in the ranges of control bam file:"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
bamToBed -i ${1} > Control.bed
mergeBed -i Control.bed > Control.merged.bed
multiBamCov -bams ${1} -bed Control.merged.bed > control_counts
awk '{ if ($4 > 0) { print } }' control_counts > filter_merged.bed
vcfintersect -b filter_merged.bed case_variants.QUAL2.filter.vcf > Case.filtered.vcf
echo "Done"
echo "Initially filtered Case-associated variants are named Case.filtered.vcf"
rm Control.bed filter_merged.bed Control.merged.bed control_counts
rm case_variants*
rm Control_initial_filter.vcf
rm *bcftools.vcf
echo ""
### Performing Somatic Variant Calling with strelka v2.9.2
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Performing Somatic Variant Calling with strelka v2.9.2:"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
echo ""
echo "for documentation, please see: https://github.com/Illumina/strelka"
echo ""
echo "downloading strelka binary from github repository"
wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
tar xvjf strelka-2.9.2.centos6_x86_64.tar.bz2
echo ""
echo "==> run demo to check successful installation"
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaSomaticWorkflowDemo.bash
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaGermlineWorkflowDemo.bash
echo "demo run was succesfull"
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Running on control and case samples: Collecting Germline variants:"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
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
printf "${CYAN}:::::::::::::::::::::\n"
echo "Variant Calling done"
end=`date +%s`
elapsed=`expr $end - $begin`
printf "${CYAN}:::::::::::::::::::::${NC}\n"
cp ./strelka_germline/results/variants/variants.vcf.gz ./strelka_germline_variants.vcf.gz
bgzip -d strelka_germline_variants.vcf.gz
grep "#" strelka_germline_variants.vcf > strelka_germline_variants_header.vcf
grep "PASS" strelka_germline_variants.vcf > strelka_germline_variants_PASS.vcf
grep -v "NoPassedVariantGTs" strelka_germline_variants_PASS.vcf > strelka_germline_variants_PASS2.vcf
cat strelka_germline_variants_header.vcf strelka_germline_variants_PASS2.vcf > strelka_germline_variants.filtered.vcf
rm strelka_germline_variants_header.vcf strelka_germline_variants_PASS.vcf strelka_germline_variants_PASS2.vcf
echo ""
echo "Fitered variants are called strelka_germline_variants.filtered.vcf"
echo ""
printf "${CYAN}:::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Continue with Somatic Variant Calling"
printf "${CYAN}:::::::::::::::::::::::::::::::::::::::::::${NC}\n"
# configuration
begin=`date +%s`
./strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${1} \
    --tumorBam ${2} \
    --referenceFasta ${3} \
    --runDir strelka_somatic
# execution on a single local machine with n parallel jobs
strelka_somatic/runWorkflow.py -m local -j ${5}
echo ""
printf "${CYAN}:::::::::::::::::::::\n"
echo "Variant Calling done"
end=`date +%s`
elapsed=`expr $end - $begin`
printf "${CYAN}:::::::::::::::::::::${NC}\n"
echo ""
printf "${CYAN}::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Filtering Germline and Somatic variants"
printf "${CYAN}::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
cp ./strelka_somatic/results/variants/somatic.snvs.vcf.gz ./strelka_somatic_variants.vcf.gz
bgzip -d strelka_somatic_variants.vcf.gz
grep "#" strelka_somatic_variants.vcf > strelka_somatic_variants_header.vcf
grep "PASS" strelka_somatic_variants.vcf > strelka_somatic_variants_PASS.vcf
cat strelka_somatic_variants_header.vcf strelka_somatic_variants_PASS.vcf > strelka_somatic_variants.filtered.vcf
rm strelka_somatic_variants_header.vcf strelka_somatic_variants_PASS.vcf
cp ./strelka_somatic/results/variants/somatic.indels.vcf.gz ./strelka_somatic_indels.vcf.gz
bgzip -d strelka_somatic_indels.vcf.gz
grep "#" strelka_somatic_indels.vcf > strelka_somatic_indels_header.vcf
grep "PASS" strelka_somatic_indels.vcf > strelka_somatic_indels_PASS.vcf
cat strelka_somatic_indels_header.vcf strelka_somatic_indels_PASS.vcf > strelka_somatic_indels.filtered.vcf
rm strelka_somatic_indels_header.vcf strelka_somatic_indels_PASS.vcf
echo ""
echo "Done"
echo ""
# Filtering Case.filtered.vcf variants file with strelka outputs
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Filtering Case.filtered.vcf variants file with strelka outputs..."
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
cat strelka_somatic_indels.filtered.vcf strelka_somatic_variants.filtered.vcf > strelka_all_somatic.vcf
grep "#" strelka_all_somatic.vcf > strelka_somatic_header.vcf
grep -v "#" strelka_all_somatic.vcf > strelka_somatic_SNVs.vcf
cat strelka_somatic_header.vcf strelka_somatic_SNVs.vcf > strelka_somatic.vcf
rm strelka_all_somatic.vcf strelka_somatic_header.vcf strelka_somatic_SNVs.vcf
vcfintersect -i strelka_germline_variants.filtered.vcf Case.filtered.vcf -r ${ref} --invert > Case.filtered.st.vcf
vcfintersect -i strelka_somatic.vcf Case.filtered.st.vcf -r ${ref} > Case.filtered.strelka.vcf
rm Case.filtered.st.vcf strelka_somatic.vcf
echo "Done"
### Annotating variants and obtaining gene list
echo ""
printf "${YELLOW}::::::::::::::::::::::::\n"
echo "==> Annotating variants"
printf "${YELLOW}::::::::::::::::::::::::${NC}\n"
bedtools intersect -a ${GTF} -b Case.filtered.strelka.vcf > Case.filtered.strelka.gtf
perl -lne 'print "@m" if @m=(/((?:gene_id)\s+\S+)/g);' Case.filtered.strelka.gtf > genes.with.variants.tabular
awk '!a[$0]++' genes.with.variants.tabular > gene_id_identifiers.tab
rm genes.with.variants.tabular
sed -i 's/gene_id //g' gene_id_identifiers.tab
sed -i 's/";//g' gene_id_identifiers.tab
sed -i 's/"//g' gene_id_identifiers.tab  
mv gene_id_identifiers.tab genes_with_variants.tabular
echo ""
### Output files
echo "All done"
echo ""
mkdir ${case_name}
rm Case.filtered.vcf
mv Case.filtered.strelka.gtf genes_with_variants.tabular Case.filtered.strelka.vcf *.filtered.vcf ./${case_name}
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "The following files are located in the the ./${case_name} folder"
echo ""
echo "(1) Case.filtered.strelka.gtf"     
echo "(2) Case.filtered.strelka.vcf"    
echo "(3) genes_with_variants.tabular"
echo "(4) strelka_germline_variants.filtered.vcf"
echo "(5) strelka_somatic_variants.filtered.vcf"
echo "(6) strelka_somatic_indels.filtered.vcf"
echo ""
echo "Corresponding to:"
echo ""
echo "(1): Case associated variants in GTF format"
echo "(2): Case-associated variants in VCF format"
echo "(3): List of genes with variants in tabular format"
echo "(4): Strelka germline variants"
echo "(5): Strelka somatic variants associated with case bam file"
echo "(6): Strelka somatic indels associated with case bam file"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"

#
} | tee logfile
#
