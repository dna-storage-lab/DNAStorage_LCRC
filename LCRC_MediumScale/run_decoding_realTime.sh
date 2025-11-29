#!/bin/bash
# set -o errexit

export PATH=bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./lib

#-------------------------

fastqpass_folder=$1

# Base available threshold
base_threshold=$2
# base_threshold=0.90

# reads num per file
target_read=4000

# Thread num for ldpc decoding
t_dec=32
# Thread number for reading files
t_file=4    # Thread files

# Threshold for index double-check (0~1), 0.8 is recommended
threshold=0.8
# 0: Correlation-based, 1: Alignment-based
mode=1

#-------------------------

result_dir=output/RealTimeDec_pool
if [ ! -d $result_dir ]; then
    mkdir -p $result_dir
fi

log_dir=$result_dir
if [ ! -d $log_dir ]; then
    mkdir -p $log_dir
fi

echo "Current time: $(date '+%Y-%m-%d %H:%M:%S')"


# Detect top folder
while true; do
    folder_name=($(find "$fastqpass_folder" -maxdepth 1 -mindepth 1 -type d | awk -F '/' '{print $NF}'))
    num_folder=${#folder_name[@]}
    if [ $num_folder -eq 2 ]; then
        chip1_folder=${folder_name[0]}
        chip2_folder=${folder_name[1]}
        echo "Experimental folder of chip 1 found: ""$fastqpass_folder"/"$chip1_folder"
        echo "Experimental folder of chip 2 found: ""$fastqpass_folder"/"$chip2_folder"
        break
    fi 
done

fastqpass_chip1_folder="$fastqpass_folder/$chip1_folder"
fastqpass_chip2_folder="$fastqpass_folder/$chip2_folder"


# Get Fastq prefix
chip1_first_fastq_suffix=_0.fastq
chip1_first_fastq_path_suffix="$fastqpass_chip1_folder"/*"$chip1_first_fastq_suffix"
while true; do
    if [ -z "$(ls -A $chip1_first_fastq_path_suffix 2>/dev/null)" ]; then
        :
    else
        chip1_first_fastq=$(ls $chip1_first_fastq_path_suffix 2>/dev/null)
        if [ -f "$chip1_first_fastq" ]; then
            chip1_fastq_prefix=$(echo "$chip1_first_fastq" | sed "s/$chip1_first_fastq_suffix\$//")
            break
        fi
    fi
done


chip2_first_fastq_suffix=_0.fastq
chip2_first_fastq_path_suffix="$fastqpass_chip2_folder"/*"$chip2_first_fastq_suffix"
while true; do
    if [ -z "$(ls -A $chip2_first_fastq_path_suffix 2>/dev/null)" ]; then
        :
    else
        chip2_first_fastq=$(ls $chip2_first_fastq_path_suffix 2>/dev/null)
        if [ -f "$chip2_first_fastq" ]; then
            chip2_fastq_prefix=$(echo "$chip2_first_fastq" | sed "s/$chip2_first_fastq_suffix\$//")
            break
        fi
    fi
done


chip1_fastq_prefix="$chip1_fastq_prefix""_"
chip2_fastq_prefix="$chip2_fastq_prefix""_"
echo "Chip-1 Fastq found:  $chip1_fastq_prefix"\*.fastq""
echo "Chip-2 Fastq found:  $chip2_fastq_prefix"\*.fastq""
echo "Current time: $(date '+%Y-%m-%d %H:%M:%S')"

###----------Identification & Error correction---------------###
Decoder_RealTime $result_dir $log_dir $chip1_fastq_prefix $chip2_fastq_prefix $t_file $t_dec $threshold $mode $base_threshold $target_read

###-------Showing images-------###
feh --geometry 270x180+1600+640 --scale-down $result_dir/Tundra.jpg &
feh --geometry 270x180+1330+640 --scale-down $result_dir/Lake.jpg &
feh --geometry 270x180+1500+850 --scale-down $result_dir/Galaxy.jpg

