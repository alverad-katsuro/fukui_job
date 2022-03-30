#!/bin/bash
#SBATCH -J 0a5b6e59
#SBATCH --time 5-00:00:00
#SBATCH -p cpu
#SBATCH --ntasks 8

date; pwd;

module load gaussian/09

echo "End job"

date

module load gaussian/09

g09 < 0a5b6e59.com > 0a5b6e59.log

date
