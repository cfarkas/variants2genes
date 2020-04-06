#!/bin/bash

control=${1}
case=${2}
ref=${3}
Rscript_path=${4}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [Control Bam File] [Case Bam File] [Reference] [path/to/Bam_Coverage.R]"
  echo ""
  echo "This program will invoke an R script to plot germline and somatic variants from Control and Case bam files"
  echo ""
  echo "[Control Bam File]: Path to control bam file"
  echo ""
  echo "[Case Bam File] [Reference]: Path to case bam file"
  echo ""
  echo "[path/to/bam_coverage.R]: path/to/bam_coverage.R, check ../scripts/bam_coverage.R"
  exit 0
fi

if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [Control Bam File] [Case Bam File] [Reference] [path/to/Bam_Coverage.R]"
  echo ""
  echo "This program will invoke an R script to plot germline and somatic variants from Control and Case bam files"
  echo ""
  echo "[Control Bam File]: Path to control bam file"
  echo ""
  echo "[Case Bam File] [Reference]: Path to case bam file"
  echo ""
  echo "[path/to/bam_coverage.R]: path/to/bam_coverage.R, check ../scripts/bam_coverage.R"
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [Control Bam File] [Case Bam File] [Reference] [path/to/Bam_Coverage.R]"
  echo ""
  echo "This program will invoke an R script to plot germline and somatic variants from Control and Case bam files"
  echo ""
  echo "[Control Bam File]: Path to control bam file"
  echo ""
  echo "[Case Bam File] [Reference]: Path to case bam file"
  echo ""
  echo "[path/to/bam_coverage.R]: path/to/bam_coverage.R, check ../scripts/bam_coverage.R"
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [Control Bam File] [Case Bam File] [Reference] [path/to/Bam_Coverage.R]"
  echo ""
  echo "This program will invoke an R script to plot germline and somatic variants from Control and Case bam files"
  echo ""
  echo "[Control Bam File]: Path to control bam file"
  echo ""
  echo "[Case Bam File] [Reference]: Path to case bam file"
  echo ""
  echo "[path/to/bam_coverage.R]: path/to/bam_coverage.R, check ../scripts/bam_coverage.R"
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [Control Bam File] [Case Bam File] [Reference] [path/to/bam_coverage.R]"; exit 1; }

if [ $# -ne 4 ]; then
  echo 1>&2 "Usage: ./`basename $0`  [Control Bam File] [Case Bam File] [Reference] [path/to/bam_coverage.R]"
  exit 3
fi
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
#
### File name definitions
control_name=$(echo "${1}" | awk -F'[.]' '{print $1}')
case_name=$(echo "${2}" | awk -F'[.]' '{print $1}')
#
### Variant Calling
echo "................................................................................."
echo ""
echo "Performing Variant Calling with SAMtools and bcftools (see: http://samtools.github.io/bcftools/):"
echo ""
echo "The output directory will be the following:"
echo ${dir1}
begin=`date +%s`
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${ref} --threads 1 -Ou ${1}| bcftools call -mv -Ov -o ${control_name}.vcf
echo "done with Control Bam file. Continue with Case bam file..."
bcftools mpileup -B -C 50 -d 250 --fasta-ref ${ref} --threads 1 -Ou ${2}| bcftools call -mv -Ov -o ${case_name}.vcf
end=`date +%s`
elapsed=`expr $end - $begin`
echo ""
echo "Variant Calling done"
echo "................................................................................."
echo Time taken: $elapsed
echo ""
### Filtering and intersecting VCF files
echo "Filtering with bcftools control and case vcf files..."
echo ""
bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3' ${control_name}.vcf > Control1.vcf
bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3' ${case_name}.vcf > Case1.vcf
vcffilter -f "QUAL > 30" Control1.vcf > Control.vcf
vcffilter -f "QUAL > 30" Case1.vcf > Case.vcf
echo ""
echo "Generating plot from Control and Case variants across chromosomes..."
echo ""
echo "Outputting graph in: ${dir1}"
echo ""
Rscript ${4}
rm ${control_name}.vcf ${case_name}.vcf Control1.vcf Case1.vcf
echo ""
echo "All Done. Check graph.pdf plot to explore variants from Control and Case bam files"
echo ""
