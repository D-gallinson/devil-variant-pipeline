#!/bin/bash

#########################################################
# Script to output basic metrics from stderr files from a
# parallel GenotypeGVCFs run.
#########################################################

path=$1

# Print header
printf "Err file\tinterval_file\tnum_intervals\tmemory\truntime\tvar/min\n"

for file in $path/*err
do
    # Get The interval file name and the number of intervals within that file
    interval_path=$(grep -m1 -o -E '\-L .*.interval_list' $file | sed 's/-L //g')
    interval_file=$(echo $interval_path | grep -o -E '/[0-9]+-scattered\.interval_list' | sed 's/\///g')
    interval_count=$(grep -c -E -v '^@' $interval_path)

    # Get total memory (reported by GenotypeGVCFs). If this is not present, set mem to 0
    mem=$(tail -n 1 $file | grep -o -P '(?<=Runtime.totalMemory\(\)=)[0-9]+')
    if [ -z "$mem" ]
    then
        mem=0
    fi
    let mem/=1024*1024

    # Get a random Variant/min near the middle of the file (very preliminary representation of the overall variant/min)
    lines=$(grep -c 'ProgressMeter' $file)
    let lines/=2
    if [[ "$lines" -eq 0 ]]
    then
        var_min=0
    else
        var_min=$(grep 'ProgressMeter' $file | tail -n +${lines} | grep -m 1 -o -E '[0-9]+\.[0-9]+$')
    fi

    # Get total runtime in hours. If this is not present, set time to 0
    time=$(grep -o -m1 -P '[0-9]*(?=\.[0-9]+ minutes)' $file)
    if [ -z "$time" ]
    then
        hour=0
        min=0
    fi
    let hour=time/60
    let min=time%60
    time="${hour}h ${min}min"

    # Print all metrics
    printf "${file}\t${interval_file}\t${interval_count}\t${mem}M\t${time}\t${var_min}\n"
done