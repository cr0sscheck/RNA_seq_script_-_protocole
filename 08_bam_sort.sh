#!/usr/bin/env bash
#SBATCH --job-name=BamSort
#SBATCH -o /data/users/amroczek/logs/BamSort_%J.out
#SBATCH -e /data/users/amroczek/logs/BamSort_%J.err
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
#sample="${directory_Data}/samplelist.tsv"
#line=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$sample")

sample_name=SRR7821939
#$(echo "$line" | awk '{print $1}')
sample_dir=Lung_WT_Cntrl
#$(echo “$line“ | awk '{print $2}')

#output directory
OUTDIR=${directory_bam_file}/${sample_dir}

#creating with samtools
apptainer exec -B /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools sort -@ 4 -m 24G -o $OUTDIR/${sample_name}_sorted.bam -T temp_bam ${directory_bam_file}/${sample_dir}/${sample_name}_mapping.bam