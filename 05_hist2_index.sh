#!/usr/bin/env bash
#SBATCH --job-name=hist2_index
#SBATCH -o logs/hist2_index.out
#SBATCH -e logs/hist2_index.err
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --mem=40GB
#SBATCH --partition=pibu_el8

#my directory
Data_dir=/data/users/amroczek/RNA_seq2
work_dir=${Data_dir}/genome_mus
unzip_genome_lland=/data/users/lland/rna_seq/ref_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa


#mkdir genome_index_mus

#unzip first 
#gunzip ${work_dir}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz


#can't use bwa because is only for fasta and we only have gtf.gz
# load sam file
apptainer exec /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2-build ${work_dir}/Mus_musculus.GRCm39.dna.primary_assembly.fa ${work_dir}/genome_index_mus

