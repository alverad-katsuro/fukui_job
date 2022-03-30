#!/bin/bash
#SBATCH -J c8bbbb7c
#SBATCH --time 5-00:00:00
#SBATCH -p cpu
#SBATCH --ntasks 8

date; pwd;

module load gaussian/09

echo "End job"

date

module load gaussian/09

g09 < c8bbbb7c.com > c8bbbb7c.log

date
