#!/usr/bin/env bash
#SBATCH --job-name=FeatCount
#SBATCH -o /data/users/amroczek/logs/FeatCount.out
#SBATCH -e /data/users/amroczek/logs/FeatCount.err
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem=4000M
#SBATCH --partition=pibu_el8


#my directory 
directory="/data/users/amroczek/RNA_seq2"
directory_Data=${directory}/first_data
directory_sam_file=${directory}/sam_datafile
directory_bam_file=${directory}/bam_datafile
#Mus_annotation_seq="/data/references/Mus_musculus/Ensembl/GRCm39/Annotation/Genes/build113/Mus_musculus.GRCm39.113.gtf"
Mus_annotation_seq=${directory}/genome_mus/Mus_musculus.GRCm39.113.gtf

#output
outdir="$directory/Count_reads_per_gene"
mkdir -p $outdir

apptainer exec -B $directory /containers/apptainer/subread_2.0.6.sif featureCounts -p -a $Mus_annotation_seq -o $outdir/count_table.txt -s 2 -T 4 $directory_bam_file/*/*_sorted.bam
