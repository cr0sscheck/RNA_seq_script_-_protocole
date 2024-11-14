#!/usr/bin/env bash
#SBATCH --job-name=hist2_try
#SBATCH --output=hist2_try.out
#SBATCH --error=hist2_try.err
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=400G
#SBATCH --partition=pibu_el8
#SBATCH --array=0-15

path1=first_data/Blood_WT_Case
path2=first_data/Blood_WT_Cntrl
path3=first_data/Lung_WT_Case
path4=first_data/Lung_WT_Cntrl

Samples=('path1/Blood_SRR7821949' 'path1/Blood_SRR7821950' 'path1/Blood_SRR7821951' 'path1/Blood_SRR7821952' 'path1/Blood_SRR7821953' 'path2/Blood_SRR7821968' 'path2/Blood_SRR7821969' 'path2/Blood_SRR7821970'
'path3/Lung_SRR7821918' 'path3/Lung_SRR7821919' 'path3/Lung_SRR7821920' 'path3/Lung_SRR7821921' 'path3/Lung_SRR7821922' 'path4/Lung_SRR7821937' 'path4/Lung_SRR7821938' 'path4/Lung_SRR7821939')
XX=${Samples[$SLURM_ARRAY_TASK_ID]}


#load modules
source /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif
source /data/users/lfalquet/SBC07107_24/scripts/module.sh

#can't use bwa because is only for fasta and we only have gtf.gz
# load sam file
exec /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2 -p 4 --dta -x genome_mus/Mus_musculus.GRCm39.113.gtf.gz -1 ${XX}_1.fastq.gz -2 ${XX}_2.fastq.gz -S ${XX}.sam