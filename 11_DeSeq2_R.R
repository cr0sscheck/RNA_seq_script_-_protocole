##############################################
#download different package
##############################################

#dowload library
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
BiocManager::install("org.Mmn.eg.db")
BiocManager::install('EnhancedVolcano')
BiocManager::install("clusterProfiler")

library("pheatmap")
library('DESeq2')
library("ggplot2")
library("EnhancedVolcano")
library('biomaRt') #to use Ensembl
library('org.Mm.eg.db')
library('clusterProfiler')

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

##############################################
#differential expression Blood
##############################################

# Extract Blood samples
blood_samples <- coldata[coldata$tissu == "Blood", ]
count_table_blood <- count_table[, blood_samples$sample]

# DDS for blood only
dds_blood <- DESeqDataSetFromMatrix(countData = count_table_blood, colData = blood_samples, design = ~ condition)
dds_blood <- DESeq(dds_blood)

# Compare Case vs Control
results_blood <- results(dds_blood, contrast = c("condition", "Case", "Control"))

# Filter and resume
DE_blood <- results_blood[which(results_blood$padj < 0.05), ]
summary(DE_blood)
head(DE_blood)

#connect ensembl mouse
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#add genes name from ensembl
ensembl <- useEnsembl(biomart = 'genes', dataset = 'mmusculus_gene_ensembl', version = 110)
ensembl_gene_ids_b <- rownames(DE_blood)

gene_info_b <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = ensembl_gene_ids_b,
                     mart = ensembl)
gene_info_b <- gene_info_b[match(rownames(DE_blood), gene_info_b$ensembl_gene_id), ]


DE_blood_with_gene_names <- cbind(gene_info_b, DE_blood)


#volcano plot to visualize DE genes
blood_volcano <- EnhancedVolcano(DE_blood_with_gene_names,
                                 lab = DE_blood_with_gene_names$external_gene_name,
                                 x = 'log2FoldChange',
                                 y = 'pvalue',
                                 col = c("pink3", "olivedrab", "cyan3", "coral1"),
                                 title = "DE genes: Blood Case vs Blood Control"
                                 , selectLab = 'Oas1a')

#explore expression of gene from the article
oas1a <-DE_blood["ENSMUSG00000052776",]
oas1a$log2FoldChange > 0 #true -- > overexpression in blood
oas1a <- plotCounts(dds_blood, "ENSMUSG00000052776", intgroup = c("condition"), returnData = TRUE)
boxplot(count ~ condition , data=oas1a, main = "Expression of Oas1a", col = c('cyan', 'green4') )


##############################################
##differential expression Lung
##############################################

# Extract Lung samples
lung_samples <- coldata[coldata$tissu == "Lung", ]
count_table_lung <- count_table[, lung_samples$sample]

# DDS for blood only
dds_lung <- DESeqDataSetFromMatrix(countData = count_table_lung, colData = lung_samples, design = ~ condition)
dds_lung <- DESeq(dds_lung)

# Compare Case vs Control
results_lung <- results(dds_lung, contrast = c("condition", "Case", "Control"))

# Filter and resume
DE_lung <- results_lung[which(results_lung$padj < 0.05), ]
summary(DE_lung)
head(DE_lung)

#connect ensembl mouse
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#add genes name from ensembl
ensembl <- useEnsembl(biomart = 'genes', dataset = 'mmusculus_gene_ensembl', version = 110)#specify version
ensembl_gene_ids_l <- rownames(DE_lung)

gene_info_l <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = ensembl_gene_ids_l,
                     mart = ensembl)
gene_info_l <- gene_info_b[match(rownames(DE_lung), gene_info_b$ensembl_gene_id), ]

DE_lung_with_gene_names <- cbind(gene_info_l, DE_lung)

#volcano plot to visualize DE genes
lung_volcano <- EnhancedVolcano(DE_lung_with_gene_names,
                                lab = DE_lung_with_gene_names$external_gene_name,
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                col = c("pink3", "olivedrab", "cyan3", "coral1"),
                                title = "DE genes: Lung Case vs Lung Control",
                                selectLab = 'Fcgr1'
)

#explore expression of gene from the article
fcgr1 <-DE_lung["ENSMUSG00000015947",]
fcgr1$log2FoldChange > 0 #True -- >overexpression in lung
fcgr1 <- plotCounts(dds_lung, "ENSMUSG00000015947", intgroup = c("condition"), returnData = TRUE)
boxplot(count ~ condition , data=fcgr1, main = "Expression of Fcgr1", col = c('cyan', 'green4') )

##############################################
# Answer question 6
##############################################
#number of gene differentially expressed in blood
num_significant_genes_inblood <- nrow(DE_blood)
print(num_significant_genes_inblood)

upregulated_blood <- DE_blood[which(DE_blood$log2FoldChange > 0), ]
downregulated_blood <- DE_blood[which(DE_blood$log2FoldChange < 0), ]

num_upregulated_blood <- nrow(upregulated_blood)
num_downregulated_blood <- nrow(downregulated_blood)

#number of gene differentially expressed in lung
num_significant_genes_inlung <- nrow(DE_lung)
print(num_significant_genes_inlung)

upregulated_lung <- DE_lung[which(DE_lung$log2FoldChange > 0), ]
downregulated_lung <- DE_lung[which(DE_lung$log2FoldChange < 0), ]

num_upregulated_lung <- nrow(upregulated_lung)
num_downregulated_lung <- nrow(downregulated_lung)

cat('the number of gene up and down regulate in lung', num_upregulated_lung, 'and', num_downregulated_lung,
      'and for the blood', num_upregulated_blood,'and', num_downregulated_blood,',respectively')
##############################################
# Overexpression Analysis
##############################################

##############################################
# Blood Case vs Blood Control 
##############################################

go_blood <- enrichGO(gene = DE_blood_with_gene_names$ensembl_gene_id, universe = names(dds), 
                     OrgDb = org.Mm.eg.db, ont= "BP", keyType = "ENSEMBL")


#dot plot: gene ratio
blood_dot <- dotplot(go_blood) + ggtitle("Blood (Case vs Control): GeneRatio")+ theme(axis.text.y = element_text(size = 9))
blood_dot
#other plots not used in the paper
#barplot: count by go terms, sorted by p-value
blood_bar <- barplot(go_blood, showCategory = 10) + ggtitle("Blood (Case vs Control): P-value")

#web plot: show relationship between go terms
goplot(go_blood, showCategory = 10) + ggtitle("GO terms - Blood Case vs Blood Control")

##############################################
# Lung Case vs Lung Control 
##############################################

go_lung <- enrichGO(gene = DE_lung_with_gene_names$ensembl_gene_id, universe = names(dds), 
                    OrgDb = org.Mm.eg.db, ont= "BP", keyType = "ENSEMBL")


#dotplot: gene ratio
lung_dot <- dotplot(go_lung ) + ggtitle("Lung (Case vs Control): GeneRatio") + theme(axis.text.y = element_text(size = 9))
lung_dot

top_terms_for_lung <- head(go_lung, n=10)  # Top 10 terms
print(as.data.frame(top_terms_for_lung))

#other plots not used in the paper
#barplot: count by go terms, sorted by p-value
lung_bar <- barplot(go_lung, showCategory = 10) + ggtitle("Lung (Case vs Control): P-value")

#web plot: show relationship between go terms
goplot(go_lung, showCategory = 10) + ggtitle("GO terms - Lung Case vs Lung Control")
