#!/usr/bin/env bash
#SBATCH --job-name=SamBam
#SBATCH -o /data/users/amroczek/logs/SamBam_%J.out
#SBATCH -e /data/users/amroczek/logs/SamBam_%J.err
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=4000M
#SBATCH --partition=pibu_el8
#SBATCH --array=1-16

#my directory 
directory="/data/users/amroczek/RNA_seq2"
directory_Data=${directory}/first_data
directory_sam_file=${directory}/sam_datafile

#data liste
sample="${directory_Data}/samplelist.tsv"
line=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$sample")

sample_name=$(echo "$line" | awk '{print $1}')
sample_dir=$(echo “$line“ | awk '{print $2}')

#out directory bam file
OUTDIR_before=${directory}/bam_datafile/${sample_dir}
mkdir -p ${OUTDIR_before}
OUTDIR=${OUTDIR_before}/${sample_name}

#creating with samtools
apptainer exec -B /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools view -hbS ${directory_sam_file}/${sample_dir}/${sample_name}_mapping.sam -o ${OUTDIR}_mapping.bam