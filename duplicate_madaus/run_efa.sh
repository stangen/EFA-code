#! /bin/bash
#
#$ -q MPI
#$ -N eccc
#$ -V
#$ -cwd
#$ -o out_files/out_$JOB_NAME$TASK_ID
#$ -e err_files/err_$JOB_NAME$TASK_ID
#$ -t 1-182

var=$JOB_NAME
#var=ecmwf #ensemble type
var2="T2M,ALT" #all variables in the netCDF, in a particular order, csv
var3="T2M,ALT" #ob types we are assimilating, csv
var4="T2M,ALT" #variable types we are updating, csv
var5=true #string, only self-update each variable = true
var6=GC #localization type, GC = gaspari cohn
var7=1000 #localization radius
var8=20130401_0000 #start date (YYYYmmdd_HHMM)

#python run_efa_script.py
python run_efa_script.py $var $var2 $var3 $var4 $var5 $var6 $var7 $var8 $SGE_TASK_ID
