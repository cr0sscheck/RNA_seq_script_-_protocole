# DE Toxoplasma - rna-seq course
In this repository, you will find all the code to reproduce my analysis based on the reads from Singhania et al. (2019) (only lung and blood tissues from wildtype and toxoplasma).

## Steps 
### Qualtiy check 
First use fastQC to assess the quality of the reads for the lung and the blood: [01_protocole_Lung](01_protocole_Lung.sh) and [02_protocole_Blood](02_protocole_Blood.sh)to create the fastqc reports, and then use MultiQC to have all your fastqc in one report with command Multiqc . , done in [03_download_gen_specie_+_multiqc](03_download_gen_specie_+_multiqc.sh).

### Creating the index
You will need to download the reference genome from ensembl [03_download_gen_specie_+_multiqc](03_download_gen_specie_+_multiqc.sh). Then run [04_hist2_index](04_hist2_index.sh), to avoid download many time, I take it in lland directory.

### Mapping the reads against the reference genome
To do that, run [05_hist2_mapping](05_hist2_mapping.sh), [06_sam_to_bam](06_sam_to_bam.sh), [07_bam_sort](07_bam_sort.sh) and [08_bam_index](08_bam_index.sh)
Then you will have a .bam and a .bam.bai file for each of your reads.

### Counting the number of reads per genes
You will need to download the annotation corresponding to the reference genome you download before. Then run [09_feature_count](09_feature_count.sh) and [10_Seq2](10_Seq2.sh) to get the count table that you will use for the next step. Download it locally.

Then you can run the R script [11_DESeq2_R](11_DESeq2_R.R).

## also find :
You have the [combineHisat2Map](combineHisat2Map.txt) to see alignment result.
