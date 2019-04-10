#!/bin/bash

control=${1}
case=${2}
threads=${3}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Threads}"
  echo ""
  echo "This script will call variants using freebayes-parallel in Control and Case bam files to obtain case-ligated variants."
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Threads}"
  echo ""
  echo "This script will call variants using freebayes-parallel in Control and Case bam files to obtain case-ligated variants"
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Threads}"
  echo ""
  echo "This script will call variants using freebayes-parallel in Control and Case bam files to obtain case-ligated variants"
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {Control Bam File} {Case Bam File} {Threads}"
  echo ""
  echo "This script will call variants using freebayes-parallel in Control and Case bam files to obtain case-ligated variants"
  echo ""
  echo "Control Bam File: File of path to Control bam file"
  echo ""
  echo "Case Bam File: File of path to Case bam file"
  echo ""
  echo "Threads: Number of CPUs for the task (integer)"
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {Control Bam file} {Case Bam file} {Threads}"; exit 1; }

if [ $# -ne 3 ]; then
  echo 1>&2 "Usage: bash ./`basename $0` {Control Bam file} {Case Bam file} {Threads}"
  exit 3
fi
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

### File name definitions
control_name=$(echo "${1}" | awk -F'[.]' '{print $1}')
case_name=$(echo "${2}" | awk -F'[.]' '{print $1}')

### Sorting bam files
echo ""
echo "Sorting Bam files ..."
samtools sort ${1} ${control_name}.sorted -@ ${threads}
samtools sort ${2} ${case_name}.sorted -@ ${threads}

echo "Done. Indexing bam files"
samtools index ${control_name}.sorted.bam
samtools index ${case_name}.sorted.bam
echo "Done"
