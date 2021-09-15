#!/bin/bash


batch=Capture1_6-11-21
logs=$WORK_BGFS/master/scripts/logs/$batch/err/9_GenotypeGVCFs/100_parallel

printf "" > parallel.times
printf "" > tmp_minutes
printf "" > tmp_interval_path
printf "" > tmp_intervals

for file in $logs
do
    