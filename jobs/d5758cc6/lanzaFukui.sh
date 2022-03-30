#!/bin/bash
#SBATCH -J d5758cc6
#SBATCH --time 5-00:00:00
#SBATCH -p cpu
#SBATCH --ntasks 8

date; pwd;

module load gaussian/09

echo "End job"

date

module load gaussian/09

g09 < d5758cc6.com > d5758cc6.log

date
