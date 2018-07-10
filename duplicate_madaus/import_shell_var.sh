#! /bin/bash

#$ -V
#$ -cwd
#$ -o out_files/practice
#$ -e err_files/practice

var='{"T2M":1,"ALT":1}'

python import_shell_var.py $var
