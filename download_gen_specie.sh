#!/usr/bin/env bash
#SBATCH --job-name=genome_mus
#SBATCH --output=genome_mus.out
#SBATCH --error=genome_mus.err
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=25G
#SBATCH --partition=pibu_el8

source /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif


mkdir /data/users/amroczek/RNA_seq2/genome_mus
cd /data/users/amroczek/RNA_seq2/genome_mus

#genome mus
wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

#annotation mus
wget ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz

#check the sum
wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/CHECKSUMS
wget ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/CHECKSUMS