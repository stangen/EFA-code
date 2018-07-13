#! /bin/bash
#
#$ -q MPI
#$ -N eccc
#$ -V
#$ -cwd
#$ -o out_files/stats_11post$JOB_NAME$TASK_ID
#$ -e err_files/stats_11post$JOB_NAME$TASK_ID
#$ -t 1-8

var=$JOB_NAME #ensemble type
var2="T2M,ALT" #all variables in the prior netCDF, in a particular order (T2M,ALT)
var3="ALT1,ALT0" #all variables in the posterior netCDF, contains ob err var
var4="ALT" #ob type that we want to run statistics on
var5="1" #observation error variance of ob type 
var6="ALT" #all obs we will have in the end, for saving the file (T2M,ALT)
var7=2013040100 #start date YYYYmmddHH
var8=2013043012 #end date YYYYmmddHH
var9=12 #hour increment between start and end date
var12=true #true for posterior, false for prior comparison
var13=1000 #localization radius
var14=none #inflation factor
var15=gridded #madis or gridded obs were assimilated?
#the 2 sge task ids are for the forecast hours to get stats on. the first one is for the starting forecast hour, 2nd is for ending forecast hour.

python mse_variance_dict_EFA.py $var $var2 $var3 $var4 $var5 $var6 $var7 $var8 $var9 $SGE_TASK_ID $SGE_TASK_ID $var12 $var13 $var14 $var15
