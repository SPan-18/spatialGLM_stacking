#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o joblog.$JOB_ID #creates a file called joblog.jobidnumber to write to.
#$ -j y
#$ -l h_data=3G,h_rt=30:00:00,highp #requests 30 minutes, 2GB of data (per core)
#$ -pe shared 1 #requests 2 cores
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system
# Notify when
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.2.2 #loads R for use
module load gcc/11.3.0
module load cmake/3.19.5

# run R code
echo 'Running h_test_mcmc.R' #prints this quote to joblog.jobidnumber
Rscript h_sptv_mcmc.R > output.$JOB_ID 2>&1
# command-line arguments: number of samples, number of repetitions, seed, command to create n samples
# outputs any text (stdout and stderr) to output.$JOB_ID