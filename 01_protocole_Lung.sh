#!/usr/bin/env bash
#SBATCH --job-name=fastqc_sample1
#SBATCH --output=fastqc_sample1.out
#SBATCH --error=fastqc_sample1.err
#SBATCH --cpus-per-task=3
#SBATCH --time=2:00:00
#SBATCH --mem=25G
#SBATCH --partition=pibu_el8
#SBATCH --array=0-2

#you have to change array between 0-2 or 0-4
#folders lung_WT_case
#folders=('SRR7821918' 'SRR7821919' 'SRR7821920' 'SRR7821921' 'SRR7821922')


#folders lung_WT_cntrl
folders=('SRR7821937' 'SRR7821938' 'SRR7821939')

XX=${folders[$SLURM_ARRAY_TASK_ID]}

#load fastqc source
source ~/containers/apptainer/fastqc-0.12.1.sif
source /data/users/lfalquet/SBC07107_24/scripts/module.sh

#create new file and go inside
mkdir /data/users/${USER}/RNA_seq2/first_data/Lung_WT_Cntrl
mkdir /data/users/${USER}/RNA_seq2/first_data/Lung_WT_Cntrl/Lung_${XX} 
cd /data/users/${USER}/RNA_seq2/first_data/Lung_WT_Cntrl/Lung_${XX}

#load the data inside my file
ln -s /data/courses/rnaseq_course/toxoplasma_de/reads/${XX}_1.fastq.gz ${XX}_1.fastq.gz
ln -s /data/courses/rnaseq_course/toxoplasma_de/reads/${XX}_2.fastq.gz ${XX}_2.fastq.gz

#fastqc quality crontrol check
fastqc -t 2 ${XX}_*.fastq.gz

#clean with fastp
fastp -i ${XX}_1.fastq.gz -I ${XX}_2.fastq.gz -o ${XX}_1trim.fastq.gz -O ${XX}_2trim.fastq.gz -j ${XX}_fastp.json -h ${XX}_fastp.html --thread 8 --trim_poly_g -l 50;