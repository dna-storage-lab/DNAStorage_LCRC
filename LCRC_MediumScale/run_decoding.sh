#!/bin/bash

export PATH=bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./lib

#-------------------------

# Fastq filename
src_fastq=$1
# src_fastq=fastq/Sample2_ONT_SUP.fastq

# Copy number of sequences
copy_num=$2

# Thread num
t_decode=$3

# Threshold for index double-check (0~1), 0.8 is recommended
threshold=$4

# 0: Correlation-based, 1: Alignment-based
mode=$5

num_refer=299700
#-------------------------

# Directory
result_dir=output/decode_pool
if [ ! -d $result_dir ]; then
    mkdir -p $result_dir
fi
log_dir=$result_dir

subfastq_prefix=$result_dir/subfastq_
# echo -e "$subfastq_prefix\n"

reads_per_file=$(awk -v a=$num_refer -v b=$copy_num -v c=$t_decode 'BEGIN{print a*b/c }')
reads_per_file=$(echo "$reads_per_file" | awk '{print int($1)}')

tot_lines=$(awk 'END{print NR}' $src_fastq)
lines_per_file=$(( reads_per_file * 4 ))
repeats=$(( tot_lines/lines_per_file/t_decode ))

echo "Total reads in raw fastq: $(( tot_lines / 4 ))"
echo "Reads per file: $reads_per_file"

################## Split FASTQ files (suffix: 05d)

split -a5 --additional-suffix=.fastq -l $lines_per_file -d $src_fastq $subfastq_prefix
echo "Raw FASTQ file ($tot_lines lines) was split into $t_decode subfiles! $(date "+%Y-%m-%d %H:%M:%S")"
echo "Lines per thread: $lines_per_file"

repeats=1
echo "Repeats: $repeats"

for(( i=0; i<repeats; ++i ));
do
    start_suffixID=$(( t_decode * i ))
    # echo "Start suffix ID $start_suffixID"

    ################## Identification & Error correction
    Decoder_OligoPool $result_dir $log_dir $subfastq_prefix $t_decode $threshold $mode $start_suffixID
done