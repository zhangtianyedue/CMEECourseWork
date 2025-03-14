#!/bin/bash
#PBS -N demographic_cluster
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-100

module load R

cd /rds/general/user/tz124/home/

Rscript tz124_HPC_2024_demographic_cluster.R

