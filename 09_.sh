#!/usr/bin/env bash
#SBATCH --job-name=SamBam
#SBATCH -o /data/users/amroczek/logs/SamBam_%J.out
#SBATCH -e /data/users/amroczek/logs/SamBam_%J.err
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=4000M
#SBATCH --partition=pibu_el8
#SBATCH --array=1-16




apptainer exec -B /containers/apptainer/subread_2.0.1–-hed695b0_0.sif