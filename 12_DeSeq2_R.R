##############################################
#download different package
##############################################
#dowload library
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2", force = T)
BiocManager::install("biomaRt")
BiocManager::install("org.Mm.eg.db")
BiocManager::install('EnhancedVolcano')
BiocManager::install("clusterProfiler")
BiocManager::install("EnhancedVolcano")

library("pheatmap")
library('DESeq2')
library("ggplot2")
library("EnhancedVolcano")
library('biomaRt')
library('matrixStats')

##############################################
#Extract data with DESeq for analyse
##############################################

# go to our table
dir_here<-getwd()
setwd('Desktop/Aurélien/Student/FriBe/1ère_master/RNA_seq/Count_reads_per_gene/')
count_table<-read.csv('table_count_analyze.txt',sep = '\t', row.names = 1, header = T) #bello dataframe

#change the name of our samples
colon_name<-c('Blood_WT_Case_1','Blood_WT_Case_2','Blood_WT_Case_3','Blood_WT_Case_4','Blood_WT_Case_5',
              'Blood_WT_Control_1', 'Blood_WT_Control_2','Blood_WT_Control_3',
              'Lung_WT_Case_1','Lung_WT_Case_2','Lung_WT_Case_3','Lung_WT_Case_4','Lung_WT_Case_5',
              'Lung_WT_Control_1','Lung_WT_Control_2','Lung_WT_Control_3')
colnames(count_table) <- colon_name

#create a df with the condition
coldata <- data.frame(
  sample = colon_name, 
  condition = c(rep("Case",5), rep("Control",3), rep("Case",5), rep("Control",3)),
  tissu= c(rep('Blood',8),rep('Lung',8))
)

#put everything in a dds data
dds <- DESeqDataSetFromMatrix(countData = count_table, colData = coldata, design = ~ tissu + condition )

#normalize data, vsd to remove the variance on the mean
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE) 

##############################################
#Plotting 
##############################################
#vizulaze the data
plotPCA(vsd, intgroup = c("condition","tissu")) + ggtitle("PCA") + theme(plot.title = element_text(hjust = 0.5))

#first visualisation
plotMA(dds)

##############################################
#try
##############################################


res_dds<-results(dds)
res_df <- as.data.frame(res_dds)
res_df$gene_id <- rownames(res_df)
#connect ensembl mouse
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Extraire les noms des protéines correspondants aux identifiants génétiques
# Supposez que vos identifiants sont dans une colonne `gene_id` du data frame `df`
protein_names <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'protein_id'),
  filters = 'ensembl_gene_id',
  values = res_df$gene_id ,
  mart = mart
)
res_df <- merge(res_df, protein_names, by.x = "gene_id", by.y = "ensembl_gene_id")

EnhancedVolcano(res_df, lab = protein_names$external_gene_name ,x='log2FoldChange',y = 'padj')

