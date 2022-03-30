#!/bin/bash
#SBATCH -J f5848555
#SBATCH --time 5-00:00:00
#SBATCH -p cpu
#SBATCH --ntasks 8

date; pwd;

module load gaussian/09

echo "End job"

date

module load gaussian/09

g09 < f5848555.com > f5848555.log

date
