#! /bin/bash
#
#$ -q MPI
#$ -N ecmwf
#$ -V
#$ -cwd
#$ -o out_files/out1000_ALT_$JOB_NAME$TASK_ID
#$ -e err_files/err1000_ALT_$JOB_NAME$TASK_ID
#$ -t 1-60

var=$JOB_NAME
#var=ecmwf #ensemble type
var2="T2M,ALT" #all variables in the netCDF, in a particular order, csv
var3="ALT,ALT" #ob types we are assimilating, csv
var4="1,0" #ob types' associated ob error variance, csv
#var3='{"ALT":0}' #ob types we are assimilating, and associated ob error varianc
var5="ALT" #variable types we are updating, csv
var6=true #string, only self-update each variable = true
var7=GC #localization type, GC = gaspari cohn
var8=1000 #localization radius
var9=20130401_0000 #start date (YYYYmmdd_HHMM)
var11="none" #inflation (>1) or 'none'
var12="gridded" #category of observations (madis or gridded)


#python run_efa_script.py
python run_efa_script.py $var $var2 $var3 $var4 $var5 $var6 $var7 $var8 $var9 $SGE_TASK_ID $var11 $var12
