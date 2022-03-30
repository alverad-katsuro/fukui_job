#!/bin/bash
#SBATCH -J 29b92efb
#SBATCH --time 5-00:00:00
#SBATCH -p cpu
#SBATCH --ntasks 8

date; pwd;

module load gaussian/09

echo "End job"

date

module load gaussian/09

g09 < 29b92efb.com > 29b92efb.log

date
