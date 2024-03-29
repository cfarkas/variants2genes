#!/bin/bash

set -e

{

usage="$(basename "$0") [-h] [-a <Control bam file>] [-b <Case bam file>] [-g <genome.fasta>] [-r <genome.gtf>] [-s <known-snps.vcf>] [-t <threads>]
This program will call variants using bcftools in Control and Case bam files to obtain case-linked variants and will filter these variants by using BEDtools and Strelka somatic variant caller.
Arguments:
    -h  show this help text
    -a  File or path to Control bam file (sorted and indexed)
    -b  File or path to Case bam file (sorted and indexed)
    -g  Reference genome (in fasta format)
    -r  Gene annotation (in GTF format)
    -s  known snps (in VCF format, for base recalibration)
    -t  Number of threads for processing (integer)"
options=':ha:b:g:r:s:t:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    a) a=$OPTARG;;
    b) b=$OPTARG;;
    g) g=$OPTARG;;
    r) r=$OPTARG;;
    s) s=$OPTARG;;
    t) t=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$a" ] || [ ! "$b" ] || [ ! "$g" ] || [ ! "$r" ] || [ ! "$s" ] || [ ! "$t" ]; then
  echo "arguments -a, -b, -g, -r, -s and -t must be provided"
  echo "$usage" >&2; exit 1
fi

# Conditions : Input existance

if [ ! -e "$a" ]; then
    echo ""
    echo "$a does not exist. Check your -a input"
    echo ""
    exit 9999 # die with error code 9999
fi

if [ ! -e "$b" ]; then
    echo ""
    echo "$b does not exist. Check your -b input"
    echo ""
    exit 9999 # die with error code 9999
fi

if [ ! -e "$g" ]; then
    echo ""
    echo "$g does not exist. Check your -g input"
    echo ""
    exit 9999 # die with error code 9999
fi

if [ ! -e "$r" ]; then
    echo ""
    echo "$r does not exist. Check your -r input"
    echo ""
    exit 9999 # die with error code 9999
fi

if [ ! -e "$s" ]; then
    echo ""
    echo "$s does not exist. Check your -s input"
    echo ""
    exit 9999 # die with error code 9999
fi

# Conditions : Getting absolute path of inputs
echo ""
a_DIR="$( cd "$( dirname "$a" )" && pwd )"
echo ""
echo "::: The absolute path of -a is $a_DIR"
echo ""
b_DIR="$( cd "$( dirname "$b" )" && pwd )"
echo ""
echo "::: The absolute path of -b is $b_DIR"
echo ""
g_DIR="$( cd "$( dirname "$g" )" && pwd )"
echo ""
echo "::: The absolute path of -g is $g_DIR"
echo ""
r_DIR="$( cd "$( dirname "$r" )" && pwd )"
echo ""
echo "::: The absolute path of -r is $r_DIR"
echo ""
s_DIR="$( cd "$( dirname "$s" )" && pwd )"
echo ""
echo "::: The absolute path of -s is $s_DIR"
echo ""

#    .---------- constant part!
#    vvvv vvvv-- the code from above
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

printf "${YELLOW}::: Defining Variables :::\n"
echo ""
echo "Defining variables:"
echo""
# https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash

FILE1="$a"
basename "$FILE1"
control_bam_file="$(basename -- $FILE1)"
control_bam_file_extension="$control_bam_file##*.}"
control_bam_file_name="${control_bam_file%.*}"
echo "The control BAM file used as input is the following: $control_bam_file"
echo "The control BAM file prefix is the following: $control_bam_file_name"
echo ""
FILE2="$b"
basename "$FILE2"
case_bam_file="$(basename -- $FILE2)"
case_bam_file_extension="$case_bam_file##*.}"
case_bam_file_name="${case_bam_file%.*}"
echo "The case BAM file used as input is the following: $case_bam_file"
echo "The case BAM file prefix is the following: $case_bam_file_name"
echo ""
FILE3="$r"
basename "$FILE3"
reference_gtf="$(basename -- $FILE3)"
echo "The reference GTF used as input is the following: $reference_gtf"
echo ""
FILE4="$g"
basename "$FILE4"
reference_genome="$(basename -- $FILE4)"
echo "The reference genome used as input is the following: $reference_genome"
echo ""
FILE5="$s"
basename "$FILE5"
known_snps="$(basename -- $FILE5)"
echo "The known SNPs collection used as input is the following: $known_snps"
echo ""
FILE6="$t"
basename "$FILE6"
threads="$(basename -- $FILE6)"
echo "The number of threads for calculation are the following: $threads"
echo ""

### Creating directory
printf "${YELLOW}::: Creating directory :::\n"
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
sec=$(date "+%Y%m%d_%H%M%S")
mkdir variants2genes_$sec
cd variants2genes_$sec
echo "Done"
echo ""

### Clearing up reference genome dict file. 
reference_genome_extension="$reference_genome##*.}"
reference_genome_name="${reference_genome%.*}"
rm -rf ${g_DIR}/${reference_genome_name}.dict

### Variant Calling
if [ ! -f ${a_DIR}/${control_bam_file}.bai ]; then
    echo " ${a}.bai file not found!. Did you forget to sort and index bam files?"
    exit 1
fi
if [ ! -f ${b_DIR}/${case_bam_file}.bai ]; then
    echo " ${b}.bai file not found!. Did you forget to sort and index bam files?"
    exit 1
fi

begin_0=`date +%s`
echo ""
printf "${YELLOW} ==== MarkDuplicates in Control and Case BAM files using picard tools ${NC}\n"
echo ""
echo ""
picard MarkDuplicates I=${a_DIR}/${control_bam_file} O=${control_bam_file_name}_marked_duplicates.bam M=${control_bam_file_name}_marked_dup_metrics.txt
picard MarkDuplicates I=${b_DIR}/${case_bam_file} O=${case_bam_file_name}_marked_duplicates.bam M=${case_bam_file_name}_marked_dup_metrics.txt
echo ""
echo "INFO ::: Done. Duplicate identification was successful. Continue."
echo ""
echo ""
printf "${YELLOW} ==== Sorting Control and Case BAM files using picard tools ${NC}\n"
echo ""
echo ""
picard SortSam I=${control_bam_file_name}_marked_duplicates.bam O=${control_bam_file_name}_marked_duplicates.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true
picard SortSam I=${case_bam_file_name}_marked_duplicates.bam O=${case_bam_file_name}_marked_duplicates.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true
echo ""
echo "INFO ::: Sorting was Done. Continue"
echo ""
echo ""
printf "${YELLOW} ==== Assigning all the reads in BAM files to a single new read-group using picard tools ${NC}\n"
echo ""
echo ""
picard AddOrReplaceReadGroups I=${control_bam_file_name}_marked_duplicates.sorted.bam O=${control_bam_file_name}.for_gatk4.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
picard AddOrReplaceReadGroups I=${case_bam_file_name}_marked_duplicates.sorted.bam O=${case_bam_file_name}.for_gatk4.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
echo ""
echo "INFO ::: Assignments were Done. Continue with gatk4 base recalibration"
echo ""
echo ""
printf "${YELLOW} ==== Performing gatk4 BaseRecalibrator on BAM files using user-provided known SNP sites in VCF format ${NC}\n"
echo ""
echo ""
gatk CreateSequenceDictionary -R ${g_DIR}/${reference_genome}
gatk IndexFeatureFile -F ${s_DIR}/${known_snps}
gatk BaseRecalibrator --input ${control_bam_file_name}.for_gatk4.bam --reference ${g_DIR}/${reference_genome} --known-sites ${s_DIR}/${known_snps} --output ${control_bam_file_name}_recal_data.table
gatk BaseRecalibrator --input ${case_bam_file_name}.for_gatk4.bam --reference ${g_DIR}/${reference_genome} --known-sites ${s_DIR}/${known_snps} --output ${case_bam_file_name}_recal_data.table
echo ""
echo "INFO ::: Recalibration was Done. We will apply base recalibration on bam files using recently created recalibration tables"
echo ""
echo ""
printf "${YELLOW} ==== Apply gatk4 BaseRecalibrator on BAM files using recalibration tables ${NC}\n"
echo ""
echo ""
gatk ApplyBQSR -R ${g_DIR}/${reference_genome} -I ${control_bam_file_name}.for_gatk4.bam --bqsr-recal-file ${control_bam_file_name}_recal_data.table -O ${control_bam_file_name}.recalibrated.bam
gatk ApplyBQSR -R ${g_DIR}/${reference_genome} -I ${case_bam_file_name}.for_gatk4.bam --bqsr-recal-file ${case_bam_file_name}_recal_data.table -O ${case_bam_file_name}.recalibrated.bam
echo ""
echo "INFO ::: Recalibration on BAM files was Done. Continue with sorting and variant calling"
echo ""
echo ""
printf "${YELLOW} ==== Sorting and index BAM files with SAMtools ${NC}\n"
echo ""
echo ""
samtools sort -o ${control_bam_file_name}.recalibrated.sorted.bam ${control_bam_file_name}.recalibrated.bam -@ ${t} && samtools index ${control_bam_file_name}.recalibrated.sorted.bam -@ ${t}
samtools sort -o ${case_bam_file_name}.recalibrated.sorted.bam ${case_bam_file_name}.recalibrated.bam -@ ${t} && samtools index ${case_bam_file_name}.recalibrated.sorted.bam -@ ${t}
echo ""
echo "INFO ::: Sorting was done, continue with variant calling"
echo ""
echo ""
printf "${YELLOW} ==== Performing Variant Calling with gatk HaplotypeCaller ${NC}\n"
echo ""
echo ""
gatk HaplotypeCaller --input ${control_bam_file_name}.recalibrated.sorted.bam --output ${control_bam_file_name}.gatk4.germline.vcf --reference ${g_DIR}/${reference_genome}
gatk HaplotypeCaller --input ${case_bam_file_name}.recalibrated.sorted.bam --output ${case_bam_file_name}.gatk4.germline.vcf --reference ${g_DIR}/${reference_genome}
echo ""
echo "INFO ::: gatk HaplotypeCaller variant calling done, continue with variant filtering"
echo ""
echo ""
printf "${YELLOW} ==== Hard-filter SNPs on multiple expressions using gatk VariantFiltration ${NC}\n"
echo ""
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
gatk VariantFiltration -V ${control_bam_file_name}.gatk4.germline.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ${control_bam_file_name}.gatk4.germline.filtered.vcf
gatk VariantFiltration -V ${case_bam_file_name}.gatk4.germline.vcf -filter "QD < 10.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ${case_bam_file_name}.gatk4.germline.filtered.vcf
# https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b
bcftools view -f PASS ${control_bam_file_name}.gatk4.germline.filtered.vcf > ${control_bam_file_name}.gatk4.germline.PASS.vcf
bcftools view -f PASS ${case_bam_file_name}.gatk4.germline.filtered.vcf > ${case_bam_file_name}.gatk4.germline.PASS.vcf
echo ""
echo "INFO ::: Filtering done, continue with bcftools variant calling"
echo ""
echo ""
printf "${YELLOW} ==== Performing Variant Calling with bcftools (see: http://samtools.github.io/bcftools/) ${NC}\n"
echo ""
echo ""
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${g_DIR}/${reference_genome} --threads ${t} -Ou ${control_bam_file_name}.recalibrated.sorted.bam | bcftools call -mv -Ov -o ${control_bam_file_name}.vcf
echo "done with Control Bam file. Continue with Case bam file..."
echo ""
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${g_DIR}/${reference_genome} --threads ${t} -Ou ${case_bam_file_name}.recalibrated.sorted.bam | bcftools call -mv -Ov -o ${case_bam_file_name}.vcf
echo ""
echo "INFO ::: bcftools variant Calling done. continue with freebayes variant calling"
echo ""
echo ""
printf "${YELLOW} ==== Performing Variant Calling with freebayes (see: https://github.com/freebayes/freebayes) ${NC}\n"
echo ""
echo ""
freebayes -f ${g_DIR}/${reference_genome} -C 3 ${control_bam_file_name}.recalibrated.sorted.bam > ${control_bam_file_name}.freebayes.vcf
echo "done with Control Bam file. Continue with Case bam file..."
echo ""
freebayes -f ${g_DIR}/${reference_genome} -C 3 ${case_bam_file_name}.recalibrated.sorted.bam > ${case_bam_file_name}.freebayes.vcf
echo ""
echo "INFO ::: Freebayes variant Calling done. Continue with variant calling filtering"
echo ""
echo ""
printf "${YELLOW} ==== Filtering and intersecting VCF files ${NC}\n"
echo ""
echo ""
### Filtering and intersecting VCF files to discover germline variants
printf "${YELLOW} ==== Filtering and intersecting VCF files to discover germline variants ${NC}\n"
echo ""
bcftools filter -e'%QUAL<10 || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 1' ${control_bam_file_name}.vcf > Control_initial_filter.vcf
bcftools filter -e'%QUAL<10 || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3' ${case_bam_file_name}.vcf > ${case_bam_file_name}.bcftools.vcf
echo "INFO ::: Done"
echo ""
echo ""
printf "${YELLOW} ==== Selecting variants in case VCF not present in control VCF ${NC}\n"
echo ""
vcfintersect -i Control_initial_filter.vcf ${case_bam_file_name}.bcftools.vcf -r ${g_DIR}/${reference_genome} --invert > case_variants.vcf
vcfintersect -i ${control_bam_file_name}.gatk4.germline.vcf ${case_bam_file_name}.gatk4.germline.PASS.vcf -r ${g_DIR}/${reference_genome} --invert > case_variants2.vcf
echo "INFO ::: Done"
echo ""
echo ""
printf "${YELLOW} ==== Case VCF file: Filtering with vcflib ---> QUAL > 20 ${NC}\n"
echo ""
vcffilter -f "QUAL > 20" case_variants.vcf > case_variants.QUAL1.filter.vcf
vcffilter -f "QUAL > 20" case_variants2.vcf > case_variants2.QUAL1.filter.vcf
echo "INFO ::: First filter done"
echo ""
echo ""
printf "${YELLOW} ==== Case VCF file: Filtering with vcflib ---> DP > 4 ${NC}\n"
echo ""
vcffilter -f "DP > 4" case_variants.QUAL1.filter.vcf > case_variants.QUAL2.filter.vcf
vcffilter -f "DP > 4" case_variants2.QUAL1.filter.vcf > case_variants2.QUAL2.filter.vcf
echo "INFO ::: Second filter done"
echo ""
echo ""
printf "${YELLOW} ==== Intersecting case variants in the ranges of control bam file ${NC}\n"
echo ""
bamToBed -i ${a_DIR}/${control_bam_file} > Control.bed
mergeBed -i Control.bed > Control.merged.bed
multiBamCov -bams ${a_DIR}/${control_bam_file} -bed Control.merged.bed > control_counts
awk '{ if ($4 > 0) { print } }' control_counts > filter_merged.bed
vcfintersect -b filter_merged.bed case_variants.QUAL2.filter.vcf > Case.filtered.vcf
vcfintersect -b filter_merged.bed case_variants2.QUAL2.filter.vcf > Case.filtered2.vcf
echo "INFO ::: Done"
echo "INFO ::: Initially filtered Case-associated variants are named Case.filtered.vcf"
rm Control.bed filter_merged.bed Control.merged.bed control_counts
rm case_variants*
rm Control_initial_filter.vcf
rm *bcftools.vcf
echo ""

# https://github.com/Illumina/strelka/releases
### Performing Somatic Variant Calling with strelka v2.9.10
echo ""
printf "${YELLOW} ==== Performing Somatic Variant Calling with strelka v2.9.10 ${NC}\n"
echo ""
echo ""
echo "INFO ::: for documentation, please see: https://github.com/Illumina/strelka"
echo ""
echo "INFO ::: downloading strelka binary from github repository"
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar xvjf strelka-2.9.10.centos6_x86_64.tar.bz2
echo ""
printf "${CYAN} ====  run demo to check successful installation\n"
bash strelka-2.9.10.centos6_x86_64/bin/runStrelkaSomaticWorkflowDemo.bash
bash strelka-2.9.10.centos6_x86_64/bin/runStrelkaGermlineWorkflowDemo.bash
printf "${CYAN} ==== INFO ::: demo run was succesfull ${NC}\n"
echo ""
echo ""
printf "${YELLOW} ==== Running on control and case samples: Collecting Germline variants ${NC}\n"
echo ""
echo ""
# configuration
begin=`date +%s`
./strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
    --bam ${a_DIR}/${control_bam_file} \
    --bam ${b_DIR}/${case_bam_file} \
    --referenceFasta ${g_DIR}/${reference_genome} \
    --runDir strelka_germline
# execution on a single local machine with n parallel jobs
strelka_germline/runWorkflow.py -m local -j ${t}
echo ""
echo "INFO ::: Variant Calling done"
echo ""
cp ./strelka_germline/results/variants/variants.vcf.gz ./strelka_germline_variants.vcf.gz
bgzip -d strelka_germline_variants.vcf.gz
grep "#" strelka_germline_variants.vcf > strelka_germline_variants_header.vcf
grep "PASS" strelka_germline_variants.vcf > strelka_germline_variants_PASS.vcf
grep -v "NoPassedVariantGTs" strelka_germline_variants_PASS.vcf > strelka_germline_variants_PASS2.vcf
cat strelka_germline_variants_header.vcf strelka_germline_variants_PASS2.vcf > strelka_germline_variants.filtered.vcf
rm strelka_germline_variants_header.vcf strelka_germline_variants_PASS.vcf strelka_germline_variants_PASS2.vcf
echo ""
echo "INFO ::: Fitered variants are called strelka_germline_variants.filtered.vcf"
echo ""
echo ""
printf "${YELLOW} ==== Continue with Somatic Variant Calling ${NC}\n"
echo ""
# configuration
begin=`date +%s`
./strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${a_DIR}/${control_bam_file} \
    --tumorBam ${b_DIR}/${case_bam_file} \
    --referenceFasta ${g_DIR}/${reference_genome} \
    --runDir strelka_somatic
# execution on a single local machine with n parallel jobs
strelka_somatic/runWorkflow.py -m local -j ${t}
echo ""
echo ""
echo "INFO ::: Variant Calling done"
echo ""
echo ""
printf "${YELLOW} ==== Filtering Germline and Somatic variants ${NC}\n"
echo ""
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
echo "INFO ::: Done"
echo ""
# Filtering Case.filtered.vcf variants file with strelka outputs
echo ""
printf "${YELLOW} ==== Filtering Case.filtered.vcf variants file with strelka outputs ${NC}\n"
echo ""
cat strelka_somatic_indels.filtered.vcf strelka_somatic_variants.filtered.vcf > strelka_all_somatic.vcf
grep "#" strelka_all_somatic.vcf > strelka_somatic_header.vcf
grep -v "#" strelka_all_somatic.vcf > strelka_somatic_SNVs.vcf
cat strelka_somatic_header.vcf strelka_somatic_SNVs.vcf > strelka_somatic.vcf
rm strelka_all_somatic.vcf strelka_somatic_header.vcf strelka_somatic_SNVs.vcf

vcfintersect -i strelka_germline_variants.filtered.vcf Case.filtered.vcf -r ${g_DIR}/${reference_genome} --invert > Case.filtered.st.vcf
vcfintersect -i strelka_germline_variants.filtered.vcf Case.filtered2.vcf -r ${g_DIR}/${reference_genome} --invert > Case.filtered2.st.vcf

vcfintersect -i strelka_somatic.vcf Case.filtered.st.vcf -r ${g_DIR}/${reference_genome} > Case.filtered.strelka.vcf
vcfintersect -i strelka_somatic.vcf Case.filtered2.st.vcf -r ${g_DIR}/${reference_genome} > Case.filtered2.strelka.vcf

echo ""
echo ""
printf "${YELLOW} ==== Combining bcftools and gatk HaplotypeCaller germline variants into a single file ${NC}\n"
echo ""
echo ""
parallel bgzip ::: Case.filtered.strelka.vcf Case.filtered2.strelka.vcf                                                             # bgzip VCF files
parallel tabix -p vcf ::: Case.filtered.strelka.vcf.gz Case.filtered2.strelka.vcf.gz                                                # tabix VCF files
bcftools merge -o bcftools-gatk-germline.vcf.gz -O z Case.filtered.strelka.vcf.gz Case.filtered2.strelka.vcf.gz  --force-samples    # join VCF files
gunzip bcftools-gatk-germline.vcf.gz
rm Case.filtered.st.vcf Case.filtered2.st.vcf strelka_somatic.vcf Case.filtered.strelka.vcf.gz Case.filtered2.strelka.vcf.gz

vcfintersect -i bcftools-gatk-germline.vcf strelka_somatic_variants.filtered.vcf -r ${g_DIR}/${reference_genome} --invert > strelka_somatic-final.vcf
vcfintersect -i bcftools-gatk-germline.vcf strelka_somatic_indels.filtered.vcf -r ${g_DIR}/${reference_genome} --invert > strelka_indels-final.vcf
echo "INFO ::: Done"
echo ""
printf "${YELLOW} ==== Joining bcftools, gatk HaplotypeCaller, freebayes and strelka2 variants to discover somatic variants with varlociraptor ${NC}\n"
echo ""
echo ""
echo "==== BGZIP, tabix and merge control and case freebayes variants:"
parallel bgzip ::: ${control_bam_file_name}.freebayes.vcf ${case_bam_file_name}.freebayes.vcf                                                                 # bgzip VCF files
parallel tabix -p vcf ::: ${control_bam_file_name}.freebayes.vcf.gz ${case_bam_file_name}.freebayes.vcf.gz                                                    # tabix VCF files
bcftools merge -o freebayes_for_varlociraptor.vcf.gz -O z ${control_bam_file_name}.freebayes.vcf.gz ${case_bam_file_name}.freebayes.vcf.gz --force-samples    # join VCF files
echo ""
echo "INFO ::: Done"
echo ""

echo "==== BGZIP, tabix and merge control and case gatk HaplotypeCaller variants:"
parallel bgzip ::: ${control_bam_file_name}.gatk4.germline.PASS.vcf ${case_bam_file_name}.gatk4.germline.PASS.vcf                                             # bgzip VCF files
parallel tabix -p vcf ::: ${control_bam_file_name}.gatk4.germline.PASS.vcf.gz ${case_bam_file_name}.gatk4.germline.PASS.vcf.gz                                # tabix VCF files
bcftools merge -o gatk_for_varlociraptor.vcf.gz -O z ${control_bam_file_name}.gatk4.germline.PASS.vcf.gz ${case_bam_file_name}.gatk4.germline.PASS.vcf.gz --force-samples    # join VCF files
echo ""
echo "INFO ::: Done"
echo ""

echo "==== BGZIP, tabix and merge control and case freebayes variants:"
parallel bgzip ::: ${control_bam_file_name}.vcf ${case_bam_file_name}.vcf                                                                    # bgzip VCF files
parallel tabix -p vcf ::: ${control_bam_file_name}.vcf.gz ${case_bam_file_name}.vcf.gz                                                       # tabix VCF files
bcftools merge -o bcftools_for_varlociraptor.vcf.gz -O z ${control_bam_file_name}.vcf.gz ${case_bam_file_name}.vcf.gz --force-samples        # join VCF files
echo ""
echo "INFO ::: Done"
echo ""

echo "==== BGZIP, tabix and merge control and case strelka2 variants:"
parallel bgzip ::: strelka_germline_variants.filtered.vcf strelka_somatic-final.vcf strelka_indels-final.vcf                                                                   # bgzip VCF files
parallel tabix -p vcf ::: strelka_germline_variants.filtered.vcf.gz strelka_somatic-final.vcf.gz strelka_indels-final.vcf.gz                                                   # tabix VCF files
bcftools merge -o strelka2_for_varlociraptor.vcf.gz -O z strelka_germline_variants.filtered.vcf.gz strelka_somatic-final.vcf.gz strelka_indels-final.vcf.gz --force-samples    # join VCF files
echo ""
echo "INFO ::: Done"
echo ""

echo "==== Removing AD tag from merged variants with bcftools"
gunzip bcftools_for_varlociraptor.vcf.gz
gunzip gatk_for_varlociraptor.vcf.gz 
gunzip freebayes_for_varlociraptor.vcf.gz 
gunzip strelka2_for_varlociraptor.vcf.gz  

bcftools annotate -x FORMAT/AD bcftools_for_varlociraptor.vcf --force > bcftools_for_varlociraptor.filter1.vcf
bcftools annotate -x FORMAT/AD freebayes_for_varlociraptor.vcf --force > freebayes_for_varlociraptor.filter1.vcf
bcftools annotate -x FORMAT/AD gatk_for_varlociraptor.vcf --force > gatk_for_varlociraptor.filter1.vcf
bcftools annotate -x FORMAT/AD strelka2_for_varlociraptor.vcf --force > strelka2_for_varlociraptor.filter1.vcf

mv bcftools_for_varlociraptor.filter1.vcf bcftools_for_varlociraptor.vcf
mv freebayes_for_varlociraptor.filter1.vcf freebayes_for_varlociraptor.vcf
mv gatk_for_varlociraptor.filter1.vcf gatk_for_varlociraptor.vcf
mv strelka2_for_varlociraptor.filter1.vcf strelka2_for_varlociraptor.vcf

printf "${YELLOW} ==== merge variants from the four variant callers ${NC}\n"
parallel bgzip ::: bcftools_for_varlociraptor.vcf freebayes_for_varlociraptor.vcf gatk_for_varlociraptor.vcf strelka2_for_varlociraptor.vcf
parallel tabix -p vcf ::: bcftools_for_varlociraptor.vcf.gz freebayes_for_varlociraptor.vcf.gz gatk_for_varlociraptor.vcf.gz strelka2_for_varlociraptor.vcf.gz   
bcftools merge -o variants_for_varlociraptor.vcf.gz -O z bcftools_for_varlociraptor.vcf.gz gatk_for_varlociraptor.vcf.gz freebayes_for_varlociraptor.vcf.gz strelka2_for_varlociraptor.vcf.gz --force-samples # join VCF files
echo ""
echo "INFO ::: Done. Continue with varlociraptor filtering"
echo ""

echo ""
printf "${YELLOW} ==== Executing varlociraptor filtering ${NC}\n"
echo ""
echo ""
printf "${CYAN} #1: Estimating alignment properties ${NC}\n"
varlociraptor estimate alignment-properties ${g_DIR}/${reference_genome} --bam ${control_bam_file_name}.recalibrated.sorted.bam > ${control_bam_file_name}.varlociraptor.alignment-properties.json
varlociraptor estimate alignment-properties ${g_DIR}/${reference_genome} --bam ${case_bam_file_name}.recalibrated.sorted.bam > ${case_bam_file_name}.varlociraptor.alignment-properties.json
echo ""
echo "INFO ::: Done"
echo ""
printf "${CYAN} #2: preprocessing variants ${NC}\n"
varlociraptor preprocess variants ${g_DIR}/${reference_genome} --alignment-properties ${control_bam_file_name}.varlociraptor.alignment-properties.json --candidates variants_for_varlociraptor.vcf.gz --bam ${control_bam_file_name}.recalibrated.sorted.bam > ${control_bam_file_name}.varlociraptor.bcf
varlociraptor preprocess variants ${g_DIR}/${reference_genome} --alignment-properties ${case_bam_file_name}.varlociraptor.alignment-properties.json --candidates variants_for_varlociraptor.vcf.gz --bam ${case_bam_file_name}.recalibrated.sorted.bam > ${case_bam_file_name}.varlociraptor.bcf
echo ""
echo ""
printf "${YELLOW} ==== Calling and filtering variants with varlociraptor ${NC}\n"
echo ""
echo ""
printf "${CYAN} #1: call variants ${NC}\n"
varlociraptor call variants tumor-normal --tumor ${case_bam_file_name}.varlociraptor.bcf --normal ${control_bam_file_name}.varlociraptor.bcf > varlociraptor-variants.bcf
bcftools convert -O v -o varlociraptor-variants.vcf varlociraptor-variants.bcf
echo ""

printf "${CYAN} #2:filter variants: somatic ${NC}\n"
varlociraptor filter-calls control-fdr varlociraptor-variants.bcf --events SOMATIC_TUMOR --fdr 0.01 --var SNV > varlociraptor-case-somatic.FDR_1e-2.bcf
bcftools convert -O v -o varlociraptor-case-somatic.FDR_1e-2.vcf varlociraptor-case-somatic.FDR_1e-2.bcf
echo ""

printf "${CYAN} #3:filter variants: germline homozygous ${NC}\n"
varlociraptor filter-calls control-fdr varlociraptor-variants.bcf --events GERMLINE_HOM --fdr 0.01 --var SNV > varlociraptor-germline_hom.FDR_1e-2.bcf
bcftools convert -O v -o varlociraptor-germline_hom.FDR_1e-2.vcf varlociraptor-germline_hom.FDR_1e-2.bcf
echo ""

printf "${CYAN} #4: filter variants: germline heterozygous ${NC}\n"
varlociraptor filter-calls control-fdr varlociraptor-variants.bcf --events GERMLINE_HET --fdr 0.01 --var SNV > varlociraptor-germline_het.FDR_1e-2.bcf
bcftools convert -O v -o varlociraptor-germline_het.FDR_1e-2.vcf varlociraptor-germline_het.FDR_1e-2.bcf

echo ""
echo "INFO ::: Filtering Done"
echo ""

echo ""
printf "${YELLOW} ==== Filtering somatic variants using final list of germline variants ${NC}\n"
echo ""

####
#### Sanity Check: Comparing combined germline variants from HaplotypeCaller and bcftools with variants from varlociraptor
####

vcfintersect -i varlociraptor-germline_het.FDR_1e-2.vcf bcftools-gatk-germline.vcf -r ${g_DIR}/${reference_genome} --invert > bcftools-gatk-germline.F1.vcf
vcfintersect -i varlociraptor-germline_hom.FDR_1e-2.vcf bcftools-gatk-germline.F1.vcf -r ${g_DIR}/${reference_genome} --invert > bcftools-gatk-germline.F2.vcf
vcfintersect -i bcftools-gatk-germline.F2.vcf varlociraptor-case-somatic.FDR_1e-2.vcf -r ${g_DIR}/${reference_genome} --invert > varlociraptor-case-somatic.FDR_1e-2.F1.vcf

mv bcftools-gatk-germline.F2.vcf case-germline.vcf
rm bcftools-gatk-germline.F1.vcf
mv varlociraptor-case-somatic.FDR_1e-2.F1.vcf varlociraptor-case-somatic.FDR_1e-2.vcf

bcftools convert -O v -o case-germline.bcf case-germline.vcf
bcftools convert -O v -o varlociraptor-case-somatic.FDR_1e-2.bcf varlociraptor-case-somatic.FDR_1e-2.vcf

echo ""
echo "INFO ::: Filtering Done"
echo ""

echo ""
printf "${YELLOW} ==== Formatting case-associated germline variants with varlociraptor specs ${NC}\n"
echo ""
#### Intersect germline variants against varlociraptor variants
# https://davetang.org/muse/2019/09/02/comparing-vcf-files/

parallel bgzip ::: case-germline.vcf varlociraptor-variants.vcf                                                            
parallel tabix -p vcf ::: case-germline.vcf.gz varlociraptor-variants.vcf.gz                                            
bcftools isec varlociraptor-variants.vcf.gz case-germline.vcf.gz -p case-germline-overlap -n =2 -w 1
cp ./case-germline-overlap/*.vcf ./case-germline-varlociraptor.vcf
rm -r -f case-germline-overlap
gunzip case-germline.vcf.gz
gunzip varlociraptor-variants.vcf.gz
echo "number of variants from bcftools output:"
grep "chr" case-germline.vcf -c
echo ""
echo "number of variants from varlociraptor intersection:"
grep "chr" case-germline-varlociraptor.vcf -c
echo ""

echo ""
echo "INFO ::: Fomatting was Done"
echo ""

echo ""
printf "${YELLOW} ==== Annotating Germline variants ${NC}\n"
echo ""
bedtools intersect -a ${r_DIR}/${reference_gtf} -b case-germline.vcf > case-germline.gtf
perl -lne 'print "@m" if @m=(/((?:gene_id)\s+\S+)/g);' case-germline.gtf > genes.with.variants.tabular
awk '!a[$0]++' genes.with.variants.tabular > gene_id_identifiers.tab
rm genes.with.variants.tabular
sed -i 's/gene_id //g' gene_id_identifiers.tab
sed -i 's/";//g' gene_id_identifiers.tab
sed -i 's/"//g' gene_id_identifiers.tab
mv gene_id_identifiers.tab genes_with_variants.tabular
echo ""

### Output files
echo ""
echo ":::: All done ::::"
echo ""

rm Case.filtered.vcf
mkdir output_files
mv case-germline.gtf genes_with_variants.tabular case-germline.vcf case-germline-varlociraptor.vcf varlociraptor-variants.vcf varlociraptor-variants.bcf varlociraptor-case-somatic.FDR_1e-2* varlociraptor-germline_h* ./output_files
printf "${CYAN}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "The following files are located in the the ./variants2genes_$sec/output_files/ folder"
echo ""
echo "(1) case-germline.gtf"
echo "(2) case-germline.vcf"
echo "(3) case-germline-varlociraptor.vcf"
echo "(4) genes_with_variants.tabular"
echo "(5) varlociraptor-variants.vcf"
echo "(6) varlociraptor-case-somatic.FDR_1e-2.vcf"
echo "(7) varlociraptor-germline_hom.FDR_1e-2.vcf"
echo "(8) varlociraptor-germline_het.FDR_1e-2.vcf"
echo ""
echo "Corresponding to:"
echo ""
echo "(1): Case-associated germline variants in GTF format (bcftools + gatk HaplotypeCaller calls)"
echo "(2): Case-associated germline variants in VCF format (bcftools + gatk HaplotypeCaller calls)"
echo "(3): Case-associated germline variants following varlociraptor specifications"
echo "(4): List of genes with germline variants in tabular format"
echo "(5): varlociraptor complete list of variants"
echo "(6): Filtered varlociraptor somatic variants (FDR<=0.01)"
echo "(7): Filtered varlociraptor common homozygous germline variants (FDR<=0.01)"
echo "(8): Filtered varlociraptor common heterozygous germline variants (FDR<=0.01)"
printf "${CYAN}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"

end_0=`date +%s`
elapsed_0=`expr $end_0 - $begin_0`
echo ""
echo Time taken: $elapsed_0
#
} | tee logfile
#
