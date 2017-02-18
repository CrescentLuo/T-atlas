#!/bin/bash
#Record current date

job_name=$1_cutAdap
echo `date` >> ${job_name}.log
echo "Run job list $job_name" >> ${job_name}.log
sleep 1

function cut_adapter {
    sleep 1
    echo `date` >> ${job_name}.log
    echo "starting job $1" >> ${job_name}.log
    fq_file_left=${1##*/}
    fq_file_right=${1##*/}
    cutadapt \
        # format of the input files \
        -f fastq \
        # Maximum number of times an adapter sequence can be removed \
        --times 2 \
        # Maximum error rate in adapter \
        -e 0.0 \
        # Minimum overlap length \
        -O 5 \
        # Minimum sequencing quality \
        --quality-cutoff 6 \
        # Minimum read length \
        -m 20 \
        # Adapters that could appear anywhere in the sequencing read (5' end or 3' end) \
        -b AAGCAGTGGTATCAACGCAGAGTACATGGG \
        -b AAGCAGTGGTATCAACGCAGAGTACT{30}VN \
        -b AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA \
        -b TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT \
        # Read 1 output \
        -o $HOME/projects/shalek2013/processed_data/S10_R1.polyATrim.adapterTrim.fastq \
        # Read 2 output \
        -p $HOME/projects/shalek2013/processed_data/S10_R2.polyATrim.adapterTrim.fastq \
        # Read 1 input \
        $HOME/projects/shalek2013/raw_data/S10_R1.fastq.gz \
        # Read 2 input \
        $HOME/projects/shalek2013/raw_data/S10_R2.fastq.gz \
        # Statistics about how many adapters were removed 
        > $HOME/projects/shalek2013/processed_data/S10_R2.polyATrim.adapterTrim.metrics
}

todo_array=($(cat $1))

max_jobs=5
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
        do_job ${todo_array[$i]} && {
            echo `date` >> ${job_name}.log
            echo finish job ${todo_array[$index]} >> ${job_name}.log
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
