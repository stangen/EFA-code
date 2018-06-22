#! /bin/bash
#
#$ -V
#$ -cwd
#$ -o tmp.out
#$ -e tmp.err
#$ -t 1-2

var=ecmwf #ensemble type
var2="T2M,ALT" #all variables in the netCDF, in a particular order
var3="T2M,ALT" #ob types we are assimilating
var4="T2M,ALT" #variable types we are updating
var5=true #string, only self-update each variable = true
var6=GC #localization type, GC = gaspari cohn
var7=1000 #localization radius

#python run_efa_script.py
python run_efa_script.py $var $var2 $var3 $var4 $var5 $var6 $var7 $SGE_TASK_ID
