#!/bin/bash
#Record current date

job_name=$1_cutAdap
if [[ -e ${job_name}.log ]];then
    rm ${job_name}.log
fi
echo `date` >> ${job_name}.log
echo "Run job list $job_name" >> ${job_name}.log
sleep 1

function cut_adapter {
    sleep 1
    echo `date` >> ${job_name}.log
    echo "starting job $1" >> ${job_name}.log
    cutadapt -f fastq --times 2 -e 0.0 -O 5 --quality-cutoff 6 -m 20 -b AAGCAGTGGTATCAACGCAGAGTACATGGG -b AAGCAGTGGTATCAACGCAGAGTACT{30}VN -b AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA -b TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -o $1_R1.polyATrim.adapterTrim.fastq -p $1_R2.polyATrim.adapterTrim.fastq $1_1.fastq.gz $1_2.fastq.gz > $1.polyATrim.adapterTrim.metrics
    sleep 5
}

todo_array=($(cat $1))

max_jobs=10
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
        cut_adapter ${todo_array[$i]} && {
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
