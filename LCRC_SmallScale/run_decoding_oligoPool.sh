#!/bin/bash
set -o errexit

export PATH=bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./lib
ccc
# Fastq filename
src_fastq=$1
# src_fastq=fastq/Sample2_ONT_SUP.fastq

# Thread num
t=$2

# Threshold for index double-check (0~1), 0.8 is recommended
threshold=$3

# 0: Correlation-based, 1: Alignment-based
mode=$4

# Thread num for ldpc decoding
t_decode=$t

work_dir=output/decode_pool
if [ ! -d $work_dir  ]; then
	mkdir -p $work_dir
fi

subfastq_prefix=$work_dir/subfastq

################## S1: Split FASTQ files (suffix: 03d) ##################
tot_lines=$(awk 'END{print NR}' $src_fastq)
lines=$(( (tot_lines+4*t-1) / (4*t) * 4 ))
split -a3 -l $lines -d $src_fastq $subfastq_prefix
echo "Raw FASTQ file ($tot_lines lines) was split into $t subfiles! $(date "+%Y-%m-%d %H:%M:%S")"
echo "Lines per thread: $lines"

consensus_file=$work_dir/consensus_cw.bin

################## S2: Identification ##################
# build bwa index
if [ $mode -eq 1 ]; then
IdentificationCons_OligoPool $work_dir $subfastq_prefix $consensus_file $t $threshold $mode 1
fi

start_time=$(date +%s%3N)
IdentificationCons_OligoPool $work_dir $subfastq_prefix $consensus_file $t $threshold $mode 0
end_time=$(date +%s%3N)
elapsed_time_identify=$((end_time - start_time))
rm $subfastq_prefix*

################## S3: Decoding ##################
ErrorCorrection_OligoPool $work_dir $consensus_file $t_decode

echo ""
echo "Runtime of read identification: $elapsed_time_identify msec"

feh --geometry 540x360+480+200 --scale-down $work_dir/img1.jpg &
feh --geometry 540x360+480+600 --scale-down $work_dir/img2.jpg
