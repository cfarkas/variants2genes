#!/bin/bash

control_bam_file=${1}
case_bam_file_=${2}
Rscript_path=${3}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {control_bam_file} {case_bam_file} {path/to/Bam_Coverage.R}"
  echo ""
  echo "This script will invoke an R script to plot coverage ranges from Control and Case bam files"
  echo ""
  echo "{control_bam_file}: Path to control bam file"
  echo ""
  echo "{case_bam_file}: Path to case bam file"
  echo ""
  echo "{path/to/bam_coverage.R}: path/to/bam_coverage.R, check ../scripts/bam_coverage.R"
  exit 0
fi

if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {control_bam_file} {case_bam_file} {path/to/Bam_Coverage.R}"
  echo ""
  echo "This script will invoke an R script to plot coverage ranges from Control and Case bam files"
  echo ""
  echo "{control_bam_file}: Path to control bam file"
  echo ""
  echo "{case_bam_file}: Path to case bam file"
  echo ""
  echo "{path/to/bam_coverage.R}: path/to/bam_coverage.R, check ../scripts/bam_coverage.R"
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {control_bam_file} {case_bam_file} {path/to/Bam_Coverage.R}"
  echo ""
  echo "This script will invoke an R script to plot coverage ranges from Control and Case bam files"
  echo ""
  echo "{control_bam_file}: Path to control bam file"
  echo ""
  echo "{case_bam_file}: Path to case bam file"
  echo ""
  echo "{path/to/bam_coverage.R}: path/to/bam_coverage.R, check ../scripts/bam_coverage.R"
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {control_bam_file} {case_bam_file} {path/to/Bam_Coverage.R}"
  echo ""
  echo "This script will invoke an R script to plot coverage ranges from Control and Case bam files"
  echo ""
  echo "{control_bam_file}: Path to control bam file"
  echo ""
  echo "{case_bam_file}: Path to case bam file"
  echo ""
  echo "{path/to/bam_coverage.R}: path/to/bam_coverage.R, check ../scripts/bam_coverage.R"
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {control_bam_file} {case_bam_file} {path/to/bam_coverage.R}"; exit 1; }

if [ $# -ne 3 ]; then
  echo 1>&2 "Usage: bash ./`basename $0`  {control_bam_file} {case_bam_file} {path/to/bam_coverage.R}"
  exit 3
fi
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
echo ""
echo "Outputting graph in: ${dir1}"
echo ""
echo "Calculating coverage for Control and Case bam files"
bamToBed -i ${1} > ${1}.bed
echo "Control file done. Continue with Case file..."
bamToBed -i ${2} > ${2}.bed
mv ${1}.bed Control_bed
mv ${2}.bed Case_bed
echo "Case file done."
echo "done"
echo "Generating plots..."
Rscript ${3}  
echo ""
echo "Done. Check graph.pdf plot"
echo ""
echo "Obtaining aligned lectures to every chromosome, haplotype and scaffolds ..."
# Obtaining aligned lectures to every chromosome, haplotype and scaffolds ...
samtools idxstats ${1} | cut -f 1,3 > ${1}.coverage
samtools idxstats ${2} | cut -f 1,3 > ${2}.coverage
join ${1}.coverage ${2}.coverage > Control_Case.coverage
sort -nrk 2,2 Control_Case.coverage > counts_per_chromosome.txt
rm ${1}.coverage ${2}.coverage Control_Case.coverage 
echo ""
echo "All Done. Check counts_per_chromosome.txt file for lecture counts in Control and Case Bam files"
