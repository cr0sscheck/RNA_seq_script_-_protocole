#!/usr/bin/env bash
#SBATCH --job-name=essay
#SBATCH -o /data/users/amroczek/logs/essay_%J.out
#SBATCH -e /data/users/amroczek/logs/essay_%J.err
#SBATCH --cpus-per-task=4
#SBATCH --time=3:00:00
#SBATCH --mem=8000M
#SBATCH --partition=pibu_el8
#SBATCH --array=1-16


#data
sample_data="/data/courses/rnaseq_course/toxoplasma_de/reads"

directory="/data/users/amroczek/RNA_seq2"
directory_Data=${directory}/first_data
work_dir="${directory}/sam_datafile"
mkdir -p "$work_dir"
Index="${directory}/genome_mus/genome_index_mus/genome_index_mus"

#data
sample="${directory_Data}/samplelist.tsv"
line=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$sample")

sample_name=$(echo "$line" | awk '{print $1}')
sample_dir=$(echo “$line“ | awk '{print $2}')
#sample_name="SRR7821921"     
#sample_dir="Lung_WT_Case"

#output directory
output_dir="${work_dir}/${sample_dir}"
mkdir -p "$output_dir"
output_file="${output_dir}/${sample_name}"

echo "mapping $sample_name in $sample_dir$" >&2

apptainer exec -B /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2 -p 4 --dta -x ${Index} -1 ${sample_data}/${sample_name}_1.fastq.gz -2 ${sample_data}/${sample_name}_2.fastq.gz -S ${output_file}_mapping.sam --rna-strandness RF 
