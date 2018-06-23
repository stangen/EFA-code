#! /bin/bash

#$ -V
#$ -cwd
#$ -o out_files/$TASK_ID
#$ -e err_files/$TASK_ID
#$ -t 10-11

var=ecmwf

python import_shell_var.py $var $SGE_TASK_ID
