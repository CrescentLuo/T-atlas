#!/bin/bash
#Record current date
task_script=$1
job_name=$2
echo $1
if [ -e ${job_name}.log ];then
    rm ${job_name}.log
fi
echo `date` "Run job list $job_name" >> ${job_name}.log
sleep 1
function Task {
    source ${task_script}
    echo `date` "starting job $1">> ${job_name}.log
}
#ls todo files first to generate job list
todo_array=($(cat $2))

max_jobs=4
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
        Task ${todo_array[$i]} && {
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
