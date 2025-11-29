#!/bin/bash
# set -o errexit

# ./compile.sh

export PATH=bin:$PATH
# export PATH=./:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./lib

#-------------------------

fastqpass_folder=$1

# Base available threshold
base_threshold=$2
# base_threshold=0.938

# reads num per file
target_read=100

# Thread num for ldpc decoding
t_dec=29
# Threshold for index double-check (0~1), 0.8 is recommended
threshold=0.8
# 0: Correlation-based, 1: Alignment-based
mode=1

#-------------------------

result_dir=output/RealTimeDec_pool
if [ ! -d $result_dir ]; then
    mkdir -p $result_dir
fi

echo "Current time: $(date '+%Y-%m-%d %H:%M:%S')"


# Get Fastq prefix
first_fastq_suffix=_0.fastq
first_fastq_path_suffix="$fastqpass_folder"/*"$first_fastq_suffix"
while true; do
    if [ -z "$(ls -A $first_fastq_path_suffix 2>/dev/null)" ]; then
        :
    else
        first_fastq=$(ls $first_fastq_path_suffix 2>/dev/null)
        if [ -f "$first_fastq" ]; then
            fastq_prefix=$(echo "$first_fastq" | sed "s/$first_fastq_suffix\$//")
            break
        fi
    fi
done

fastq_prefix="$fastq_prefix""_"
echo "Fastq found:  $fastq_prefix"\*.fastq""
echo "Current time: $(date '+%Y-%m-%d %H:%M:%S')"

###----------Identification & Error correction---------------###
Decoder_RealTime $result_dir $fastq_prefix $t_dec $threshold $mode $base_threshold $target_read

###-------Showing images-------###
feh --borderless --geometry 540x360+900+640 --scale-down $result_dir/img1.jpg &
feh --borderless --geometry 540x360+1350+690 --scale-down $result_dir/img2.jpg
