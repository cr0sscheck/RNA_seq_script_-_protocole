#!/usr/bin/env bash
#SBATCH --job-name=CountTable
#SBATCH -o /data/users/amroczek/logs/CountTable.out
#SBATCH -e /data/users/amroczek/logs/CountTable.err
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=4000M
#SBATCH --partition=pibu_el8

#my directory 
directory="/data/users/amroczek/RNA_seq2"
directory_bam_file=${directory}/bam_datafile
count_table=${directory}/Count_reads_per_gene/count_table.txt
outdir="$directory/Count_reads_per_gene"
mkdir -p $outdir
table_count_analyse="table_count_analyze.txt"

# define an array of sample names and corresponding file paths
declare -a SAMPLES=(
    "SRR7821918"
    "SRR7821919"
    "SRR7821920"
    "SRR7821921"
    "SRR7821922"
    "SRR7821937"
    "SRR7821938"
    "SRR7821939"
    "SRR7821949"
    "SRR7821950"
    "SRR7821951"
    "SRR7821952"
    "SRR7821953"
    "SRR7821968"
    "SRR7821969"
    "SRR7821970"
)

#extract relevant column
awk -F'\t' 'NR > 1 {
    printf "%s\t", $1; 
    for (i = 7; i <= NF; i++) {
        printf "%s", $i; 
        if (i < NF) printf "\t";
    }
    printf "\n";
}' "${count_table}" > "${outdir}/$table_count_analyse"

for SAMPLE in "${SAMPLES[@]}"; do
    sed -i "s|${directory_bam_file}/*/${SAMPLE}sorted.bam|${SAMPLE}|" "${outdir}/$table_count_analyse"
done