#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000MB
#SBATCH --time=01:30:00
#SBATCH --job-name=reads_per_gene
#SBATCH -o /data/users/amroczek/logs/read_per_genes_%j.out
#SBATCH -e /data/users/amroczek/logs/read_per_genes_%j.err
#SBATCH --partition=pibu_el8

#my directory 
directory="/data/users/amroczek/RNA_seq2"
directory_Data=${directory}/first_data
Output_directory=${directory}/Count_reads_per_gene

touch "$OUTPUTDIR"/number_read.txt


## Loop through all .fastq files in the directory
for file in ${directory_Data}/*/*/*.fastq.gz; do
    # Count the number of 4-line groups in the file
    read_count=$(zcat "$file" | wc -l)
    read_count=$((read_count / 4))

    # Append the result to the output file
    echo "$file has $read_count reads" >> "$OUTPUTDIR"/number_read.txt
done
