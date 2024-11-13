#!/usr/bin/env bash
#SBATCH --job-name=genome_mus
#SBATCH --output=genome_mus.out
#SBATCH --error=genome_mus.err
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=25G
#SBATCH --partition=pibu_el8


#hist2 tools sam tools