#!/bin/bash
#SBATCH -J circRNAs
#SBATCH -o output.out
#SBATCH -e error.out
#SBATCH -t 5:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

module load bioinfo/bedtools2-2.29.0

# We make sure here that the circrna env is active
#source circrnaenv/bin/activate

#snakemake -p --jobs 8 -n
snakemake -p --jobs 8
