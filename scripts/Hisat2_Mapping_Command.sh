#!/bin/bash
#Record current date

job_name=$1
if [[ -e ${job_name}.log]]{
    rm ${job_name}.log
}
echo `date` "Run job list $job_name" >> ${job_name}.log
sleep 1
function tophat2_mapping_paired_end {
    echo `date` >> ${job_name}.log
    echo "starting job $1" >> ${job_name}.log
    fastq_file=${1##*/}
    tophat2 -g 1 -a 6 --microexon-search -p 10 -G /picb/rnomics1/database/Mouse/mm10/annotation/knownGene.gtf -o ./${fastq_file%%.*}  /picb/rnomics1/database/Mouse/mm10/genome/mm10_all $1 &> ${fastq_file%%.*}_tophat2.log
    sleep 1
}
function hisat2_mapping_paired_end {
    echo `date` "starting job $1" >> ${job_name}.log
    fastq_file=${1##*/}
    fastq_file=${fastq_file%%.*}
    # -x hisat2 index dir
    # -1 -2 paired-end fastq files
    touch ${fastq_file}_ht2.log
    echo $fastq_file
    hisat2 -x /picb/rnomics1/database/Mouse/mm10/genome/mm10_all_index --known-splicesite-infile /picb/rnomics1/database/Mouse/mm10/annotation/ref_all_spsites.txt -1 $1_R1.polyATrim.adapterTrim.fastq -2 $1_R2.polyATrim.adapterTrim.fastq -S ${fastq_file}_ht2.sam -p 10 &> ${fastq_file}_ht2.log
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
        hisat2_mapping_paired_end ${todo_array[$i]} && {
            echo `date` "finish job ${todo_array[$i]}" >> ${job_name}.log
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
