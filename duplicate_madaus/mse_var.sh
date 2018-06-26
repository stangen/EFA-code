#! /bin/bash
#
#$ -q reg
#$ -l h=enkf8
#$ -N eccc
#$ -V
#$ -cwd
#$ -o out_files/stats_$JOB_NAME$TASK_ID
#$ -e err_files/stats_$JOB_NAME$TASK_ID
#$ -t 1-8

var=$JOB_NAME #ensemble type
var2="T2M,ALT" #all variables in the netCDF, in a particular order
var3="ALT" #ob type that we want to run statistics on
var4=2013040100 #start date YYYYmmddHH
var5=2013063012 #end date YYYYmmddHH
var6=12 #hour increment between start and end date
var9=false #true for posterior, false for prior comparison
#the 2 sge task ids are for the forecast hours to get stats on. the first one is for the starting forecast hour, 2nd is for ending forecast hour.

python mse_variance_dict_EFA.py $var $var2 $var3 $var4 $var5 $var6 $SGE_TASK_ID $SGE_TASK_ID $var9
