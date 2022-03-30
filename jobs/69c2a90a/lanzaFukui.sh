#!/bin/bash
#SBATCH -J 69c2a90a
#SBATCH --time 5-00:00:00
#SBATCH -p cpu
#SBATCH --ntasks 8

date; pwd;

module load gaussian/09

echo "End job"

date

module load gaussian/09

g09 < 69c2a90a.com > 69c2a90a.log

date
