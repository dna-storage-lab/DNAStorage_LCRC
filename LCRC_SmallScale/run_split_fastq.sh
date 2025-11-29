#!/bin/bash

# fastq file
fastq_file=$1
# thread number
t=$2
# sub. fastq prefix
subfastq_prefix=$3

tot_lines=$(awk 'END{print NR}' $fastq_file)
lines=$(( (tot_lines+4*t-1) / (4*t) * 4 ))
# Split FASTQ files (suffix: 03d)
split -a3 -l $lines -d $fastq_file $subfastq_prefix

echo "Raw FASTQ file ($tot_lines lines) was split into $t subfiles! $(date "+%Y-%m-%d %H:%M:%S")"
echo "Lines per thread: $lines"
