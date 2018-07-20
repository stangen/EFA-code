#! /bin/bash
#
#$ -q MPI
# -l h=enkf9
#$ -N eccc
#$ -V
#$ -cwd
#$ -o out_files/stats_100$JOB_NAME$TASK_ID
#$ -e err_files/stats_100$JOB_NAME$TASK_ID
#$ -t 1-8

var=$JOB_NAME #ensemble type
var2="T2M,ALT" #all variables in the prior netCDF, in a particular order (T2M,ALT)
var3="ALT1,ALT0,ALT0-1" #all variables in the posterior netCDF, contains ob err var
var4="ALT" #ob type that we want to run statistics on
var5="0.1" #observation error variance of ob type 
var6="ALT" #all obs we will have in the end, for saving the file (T2M,ALT)
var7=2013040100 #start date YYYYmmddHH
var8=2013043012 #end date YYYYmmddHH
var9=12 #hour increment between start and end date
var12=true #true for posterior, false for prior comparison
var13=100 #localization radius
var14=none #inflation factor
var15=gridded #madis or gridded obs were assimilated?
var16="true" #new format: true = new filename saving conventions, false= no
var17="54" #end forecast hour of forecast
var18="-180,180,90,0,3" #left, right, top, bottom, spacing grid dimensions
var19="true" #use observation error variance in netCDF filename/var names
var20="False" #were observations assimilated to self-update their corresponding variables? 
#the 2 sge task ids are for the forecast hours to get stats on. the first one is for the starting forecast hour, 2nd is for ending forecast hour.

python mse_variance_gridded.py $var $var2 $var3 $var4 $var5 $var6 $var7 $var8 $var9 $SGE_TASK_ID $SGE_TASK_ID $var12 $var13 $var14 $var15 $var16 $var17 $var18 $var19 $var20
