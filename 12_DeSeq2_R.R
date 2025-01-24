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
library('matrixStats')#(?)
library('org.Mm.eg.db')

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
oas1a <- plotCounts(dds, "ENSMUSG00000052776", intgroup = c("condition"), returnData = TRUE)
boxplot(count ~ condition , data=oas1a, main = "Expression of Oas1a", col = c('cyan', 'green4') )

##############################################
##differential expression all
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

#extract normalize value 
mat <- assay(vsd)
rownames(mat) <- ifelse(
  rownames(mat) %in% protein_names$ensembl_gene_id, 
  protein_names$external_gene_name[match(rownames(mat), protein_names$ensembl_gene_id)], 
  rownames(mat)
)
#annotate colonne
annotation_col<- coldata[, c("tissu", "condition")]
rownames(annotation_col) <- coldata$sample


top_var_genes <- head(order(rowVars(mat), decreasing = TRUE), 50) 
mat_subset <- mat[top_var_genes, ]
summary(mat_subset)
head(mat_subset)
pheatmap(mat_subset,
         annotation_col = annotation_col ,
         cluster_rows = TRUE,            
         cluster_cols = TRUE,            
         scale = "row",                 
         fontsize_row = 8,         
         fontsize_col = 10
)

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
ensembl <- useEnsembl(biomart = 'genes', dataset = 'mmusculus_gene_ensembl', version = 110)
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
fcgr1 <- plotCounts(dds, "ENSMUSG00000015947", intgroup = c("condition"), returnData = TRUE)
boxplot(count ~ condition , data=fcgr1, main = "Expression of Fcgr1", col = c('cyan', 'green4') )

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

##############################################
# Lung Case vs Lung Control 
##############################################

go_lung <- enrichGO(gene = DE_lung_with_gene_names$ensembl_gene_id, universe = names(dds), 
                    OrgDb = org.Mm.eg.db, ont= "BP", keyType = "ENSEMBL")


#dotplot: gene ratio
lung_dot <- dotplot(go_lung ) + ggtitle("Lung (Case vs Control): GeneRatio") + theme(axis.text.y = element_text(size = 9))
lung_dot
