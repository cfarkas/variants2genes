#!/bin/bash

set -e 

usage="$(basename "$0") [-h] [-a <Control bam file>] [-b <Case bam file>] [-g <genome.fasta>] [-p <path/to/Bam_Coverage.R>]
This program will invoke an R script to plot filtered variants from Control and Case bam files, for genome-wide inspection.
Arguments:
    -h  show this help text
    -a  File or path to Control bam file (sorted and indexed)
    -b  File or path to Case bam file (sorted and indexed)
    -g  Reference genome (in fasta format)
    -p  path/to/bam_coverage.R. Check varianst2genes/scripts/bam_coverage.R"
options=':ha:b:g:p:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    a) a=$OPTARG;;
    b) b=$OPTARG;;
    g) g=$OPTARG;;
    p) p=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$a" ] || [ ! "$b" ] || [ ! "$g" ] || [ ! "$p" ]; then
  echo "arguments -a, -b, -g and -p must be provided"
  echo "$usage" >&2; exit 1
fi

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

#    .---------- constant part!
#    vvvv vvvv-- the code from above
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

### File name definitions
control_name=$(echo "${a}" | awk -F'[.]' '{print $1}')
case_name=$(echo "${b}" | awk -F'[.]' '{print $1}')

### Variant Calling
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "==> Performing Variant Calling with bcftools (see: http://samtools.github.io/bcftools/):"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
echo ""
begin=`date +%s`
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${g} --threads 1 -Ou ${a}| bcftools call -mv -Ov -o ${control_name}.vcf
echo "done with Control Bam file. Continue with Case bam file..."
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${g} --threads 1 -Ou ${b}| bcftools call -mv -Ov -o ${case_name}.vcf
end=`date +%s`
elapsed=`expr $end - $begin`
echo ""
printf "${CYAN}:::::::::::::::::::::\n"
echo "Variant Calling done"
echo ""
echo Time taken: $elapsed
printf "${CYAN}:::::::::::::::::::::${NC}\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo "Filtering with bcftools control and case vcf files..."
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
echo ""
bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3' ${control_name}.vcf > Controlv.vcf
bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3' ${case_name}.vcf > Casev.vcf
vcffilter -f "QUAL > 30" Controlv.vcf > Control.vcf
vcffilter -f "QUAL > 30" Casev.vcf > Case.vcf
echo ""
echo "Generating plot from Control and Case variants across chromosomes..."
echo ""
echo "Outputting graph in: ${dir1}"
echo ""
Rscript ${p}
rm ${control_name}.vcf ${case_name}.vcf Controlv.vcf Casev.vcf
printf "${CYAN}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
echo ""
echo "All Done. Check graph.pdf plot to explore variants from Control and Case bam files"
echo ""
printf "${CYAN}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
