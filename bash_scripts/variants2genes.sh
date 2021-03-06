#!/bin/bash

set -e

{

usage="$(basename "$0") [-h] [-a <Control bam file>] [-b <Case bam file>] [-g <genome.fasta>] [-r <genome.gtf>] [-t <threads>]
This program will call variants using bcftools in Control and Case bam files to obtain case-linked variants and will filter these variants by using BEDtools and Strelka somatic variant caller.
Arguments:
    -h  show this help text
    -a  File or path to Control bam file (sorted and indexed)
    -b  File or path to Case bam file (sorted and indexed)
    -g  Reference genome (in fasta format)
    -r  Gene annotation (in GTF format)
    -t  Number of threads for processing (integer)"
options=':ha:b:g:r:t:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    a) a=$OPTARG;;
    b) b=$OPTARG;;
    g) g=$OPTARG;;
    r) r=$OPTARG;;
    t) t=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$a" ] || [ ! "$b" ] || [ ! "$g" ] || [ ! "$r" ] || [ ! "$t" ]; then
  echo "arguments -a, -b, -g, -r and -t must be provided"
  echo "$usage" >&2; exit 1
fi

#    .---------- constant part!
#    vvvv vvvv-- the code from above
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
begin=`date +%s`
echo "Cleaning directory..."
rm -rf strelka*
echo ""
### File name definitions
control_name=$(echo "${a}" | awk -F'[.]' '{print $1}')
case_name=$(echo "${b}" | awk -F'[.]' '{print $1}')
### Variant Calling
if [ ! -f ${a}.bai ]; then
    echo " ${a}.bai file not found!. Did you forget to sort and index bam files?"
    exit 1
fi
if [ ! -f ${b}.bai ]; then
    echo " ${b}.bai file not found!. Did you forget to sort and index bam files?"
    exit 1
fi
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Performing Variant Calling with bcftools (see: http://samtools.github.io/bcftools/):"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
echo ""
begin=`date +%s`
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${g} --threads ${t} -Ou ${a}| bcftools call -mv -Ov -o ${control_name}.vcf
echo "done with Control Bam file. Continue with Case bam file..."
echo ""
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${g} --threads ${t} -Ou ${b}| bcftools call -mv -Ov -o ${case_name}.vcf
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
vcfintersect -i Control_initial_filter.vcf ${case_name}.bcftools.vcf -r ${g} --invert > case_variants.vcf
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
bamToBed -i ${a} > Control.bed
mergeBed -i Control.bed > Control.merged.bed
multiBamCov -bams ${a} -bed Control.merged.bed > control_counts
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
    --bam ${a} \
    --bam ${b} \
    --referenceFasta ${g} \
    --runDir strelka_germline
# execution on a single local machine with n parallel jobs
strelka_germline/runWorkflow.py -m local -j ${t}
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
    --normalBam ${a} \
    --tumorBam ${b} \
    --referenceFasta ${g} \
    --runDir strelka_somatic
# execution on a single local machine with n parallel jobs
strelka_somatic/runWorkflow.py -m local -j ${t}
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
vcfintersect -i strelka_germline_variants.filtered.vcf Case.filtered.vcf -r ${g} --invert > Case.filtered.st.vcf
vcfintersect -i strelka_somatic.vcf Case.filtered.st.vcf -r ${g} > Case.filtered.strelka.vcf
rm Case.filtered.st.vcf strelka_somatic.vcf
vcfintersect -i Case.filtered.strelka.vcf strelka_somatic_variants.filtered.vcf -r ${g} --invert > somatic-final.vcf
vcfintersect -i Case.filtered.strelka.vcf strelka_somatic_indels.filtered.vcf -r ${g} --invert > indels-final.vcf
echo "Done"
### Annotating variants and obtaining gene list
echo ""
printf "${YELLOW}::::::::::::::::::::::::\n"
echo "==> Annotating variants"
printf "${YELLOW}::::::::::::::::::::::::${NC}\n"
bedtools intersect -a ${r} -b Case.filtered.strelka.vcf > Case.filtered.strelka.gtf
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
mv Case.filtered.strelka.gtf genes_with_variants.tabular Case.filtered.strelka.vcf *filtered.vcf somatic-final.vcf indels-final.vcf ./${case_name}
printf "${CYAN}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "The following files are located in the the ./${case_name} folder"
echo ""
echo "(1) Case.filtered.strelka.gtf"
echo "(2) Case.filtered.strelka.vcf"
echo "(3) genes_with_variants.tabular"
echo "(4) strelka_germline_variants.filtered.vcf"
echo "(5) strelka_somatic_variants.filtered.vcf"
echo "(6) strelka_somatic_indels.filtered.vcf"
echo "(7) somatic-final.vcf"
echo "(8) indels-final.vcf"
echo ""
echo "Corresponding to:"
echo ""
echo "(1): Case associated variants in GTF format"
echo "(2): Case-associated variants in VCF format"
echo "(3): List of genes with variants in tabular format"
echo "(4): Strelka germline variants"
echo "(5): Strelka somatic variants associated with case bam file"
echo "(6): Strelka somatic indels associated with case bam file"
echo "(7): Filtered somatic variants associated with case bam file"
echo "(8): Filtered indels variants associated with case bam file"
printf "${CYAN}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"

#
} | tee logfile
#
