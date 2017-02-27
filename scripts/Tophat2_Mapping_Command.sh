#!/bin/bash
#Record current date

job_name=$1
echo `date` >> ${job_name}.log
echo "Run job list $job_name" >> ${job_name}.log
sleep 1
function tophat2_mapping_paired_end {
    echo `date` >> ${job_name}.log
    echo "starting job $1" >> ${job_name}.log
    fastq_file=${1##*/}
    tophat2 -g 1 -a 6 --microexon-search -p 10 -G /picb/rnomics1/database/Mouse/mm10/annotation/knownGene.gtf -o ./${fastq_file%%.*}  /picb/rnomics1/database/Mouse/mm10/genome/mm10_all $1 &> ${fastq_file%%.*}_tophat2.log
    sleep 1
}
function hisat2_mapping_paired_end {
    echo `date` >> ${job_name}.log
    echo "starting job $1" >> ${job_name}.log
    fastq_file=${1}
    # -x hisat2 index dir
    # -1 -2 paired-end fastq files
    hisat2 -x hg19_index -1 light-1_R1.fq -2 light-1_R2.fq -S light-1.sam -p 16 2>light-1.log
}
#ls todo files first to generate job list
todo_array=($(cat $1))

max_jobs=2
FIFO_FILE="./$$.fifo"
mkfifo $FIFO_FILE
exec 6<>$FIFO_FILE
rm -f $FIFO_FILE

for((i=1; i<=$max_jobs; i++));do
    echo
done >&6

for((i=0; i<${#todo_array[*]}; i++));do
    read -u6
    {
        tophat2_mapping_paired_end ${todo_array[$i]} && {
            echo `date` >> ${job_name}.log
            echo finish job ${todo_array[$i]} >> ${job_name}.log
            sleep 1
        } || {
            echo job ${todo_array[$i]} error >> ${job_name}.log
        }
        echo >&6
    }& 
done

# wati for all jobs to complete
wait
exec 6>&-
