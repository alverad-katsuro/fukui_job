#!/bin/bash
#SBATCH -J a9628903
#SBATCH --time 5-00:00:00
#SBATCH -p cpu
#SBATCH --ntasks 8

date; pwd;

module load gaussian/09

echo "End job"

date

module load gaussian/09

g09 < a9628903.com > a9628903.log

date
