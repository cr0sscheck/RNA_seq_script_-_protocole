#!/usr/bin/env bash
#SBATCH --job-name=BamIndex
#SBATCH -o /data/users/amroczek/logs/BamIndex_%J.out
#SBATCH -e /data/users/amroczek/logs/BamIndex_%J.err
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=30000M
#SBATCH --partition=pibu_el8
##SBATCH --array=1-16

#my directory 
directory="/data/users/amroczek/RNA_seq2"
directory_Data=${directory}/first_data
directory_sam_file=${directory}/sam_datafile
directory_bam_file=${directory}/bam_datafile

#data
sample="${directory_Data}/samplelist.tsv"
line=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$sample")

sample_name=SRR7821939
sample_dir=Lung_WT_Cntrl
#sample_name=$(echo "$line" | awk '{print $1}')
#sample_dir=$(echo “$line“ | awk '{print $2}')

#output directory
OUTDIR=${directory_bam_file}/${sample_dir}/${sample_name}

#creating with samtools
apptainer exec -B /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools index ${OUTDIR}_sorted.bam