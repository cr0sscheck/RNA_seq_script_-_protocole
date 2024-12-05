#dowload library
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2", force = T)
BiocManager::install("biomaRt")
BiocManager::install("org.Mm.eg.db")
BiocManager::install('EnhancedVolcano')
BiocManager::install("clusterProfiler")

library("pheatmap")
library('DESeq2')
library("tximport")
library("readr")
library("tximportData")
library("ggplot2")

# go to our table
dir_here<-getwd()
setwd('Desktop/Aurélien/Student/FriBe/1ère_master/RNA_seq/Count_reads_per_gene/')
count_table<-read.csv('table_count_analyze.txt',sep = '\t', row.names = 1, header = T) #bello dataframe

colon_name<-c('Blood_WT_Case_1','Blood_WT_Case_2','Blood_WT_Case_3','Blood_WT_Case_4','Blood_WT_Case_5',
              'Blood_WT_Control_1', 'Blood_WT_Control_2','Blood_WT_Control_3',
              'Lung_WT_Case_1','Lung_WT_Case_2','Lung_WT_Case_3','Lung_WT_Case_4','Lung_WT_Case_5',
              'Lung_WT_Control_1','Lung_WT_Control_2','Lung_WT_Control_3')

coldata <- data.frame(
  sample = colon_name, 
  condition = c(rep("Case",5), rep("Control",3), rep("Case",5), rep("Control",3)),
  tissu= c(rep('Blood',8),rep('Lung',8))
)

colnames(count_table) <- coldata$sample
coldata$tissu <- factor(coldata$tissu)

dds <- DESeqDataSetFromMatrix(countData = count_table, colData = coldata, design = ~ tissu + condition )

#normalize data
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = c("condition","tissu")) + ggtitle("PCA") + theme(plot.title = element_text(hjust = 0.5))

head(results(dds))
rlogged_dds<-rlog(dds)
normalized_counts <- counts(dds, normalized = TRUE)  

#first visualisation
plotMA(dds)

res_dds<-results(dds)
