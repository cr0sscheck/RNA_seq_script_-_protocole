#!/usr/bin/env bash
#SBATCH --job-name=hist2_index
#SBATCH --output=hist2_index.out
#SBATCH --error=hist2_index.err
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=40G
#SBATCH --partition=pibu_el8



#load modules
#source /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif
source /data/users/lfalquet/SBC07107_24/scripts/module.sh

#can't use bwa because is only for fasta and we only have gtf.gz
# load sam file
apptainer exec /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2-build genome_mus/Mus_musculus.GRCm39.113.gtf.gz genome_mus/mus_genome_index