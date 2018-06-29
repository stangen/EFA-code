#! /bin/bash
#
#$ -q MPI
#$ -N ecmwf
#$ -V
#$ -cwd
#$ -o out_files/stats_loc500$JOB_NAME$TASK_ID
#$ -e err_files/stats_loc500$JOB_NAME$TASK_ID
#$ -t 1-8

var=$JOB_NAME #ensemble type
var2="T2M,ALT" #all variables in the netCDF, in a particular order
var3="ALT" #ob type that we want to run statistics on
var4="ALT" #all obs we will have in the end, for saving the file (T2M,ALT)
var5=2013040100 #start date YYYYmmddHH
var6=2013043012 #end date YYYYmmddHH
var7=12 #hour increment between start and end date
var10=true #true for posterior, false for prior comparison
var11=1000
#the 2 sge task ids are for the forecast hours to get stats on. the first one is for the starting forecast hour, 2nd is for ending forecast hour.

python mse_variance_dict_EFA.py $var $var2 $var3 $var4 $var5 $var6 $var7 $SGE_TASK_ID $SGE_TASK_ID $var10 $var11
