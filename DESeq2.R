# Gene-level differential expression analysis using DESeq2, including RIN as covariate
##Load libraries
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(GenomicFeatures)
library(biomaRt)
library(ggplot2)
library(scales)
library(viridis)
library(apeglm)

#Load count table (merged count table from featurecount output)
countdata <- read.table("ccRCC_PCRcDNA_counts.txt", header = TRUE, sep = '\t', row.names = 1)
head(countdata)

#turn table into matrix
countdata <- as.matrix(countdata)

#Load metadata - columns: sample, condition - see example from github folder
samps <- read.table("samples.txt", header = TRUE, sep = '\t', row.names = 1)
head(samps)

#Set up factors
samps$condition <- factor(samps$condition)


#set up one dds without covariate
dds <- DESeqDataSetFromMatrix(countData=round(countdata), colData=samps, design=~condition)
dds

#Set up factor levels
dds$condition <- factor(dds$condition, levels = c("nonrecurrent","recurrent"))

#prefiltering - keep only rows that have at least 1 reads across samples
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

#Run DESeq2
dds <- DESeq(dds)

#total number of normalised count per sample - should be similar across samples
colSums(counts(dds, normalized=T))

#plot dispersion plots
plotDispEsts(dds, main="PCS_Dispersion_plot")

#transform dds results to generate pca plot
rld <- rlogTransformation(dds)

#create PCA plots
DESeq2::plotPCA(rld, intgroup="condition")

#Calculate DESeq2 results
res <- results(dds)

#produce simple MA plots
plotMA(res)

#More advanced MA plot - padj<.1 - ggplot
res.MA <- as.data.frame(res)
res.MA$significant <- ifelse(res.MA$padj < .1, "Significant", NA)
ggplot(res.MA,aes(baseMean, log2FoldChange, colour=padj))+ geom_point(size=1)+ scale_y_continuous(limits=c(-3, 3), oob=scales::squish)+ scale_x_log10()+geom_hline(yintercept = 0, colour="darkorchid4", linewidth=1, linetype="longdash")+ labs(x="mean of normalized counts", y="log fold change")+ scale_colour_viridis(direction=-1, trans='sqrt')+ theme_bw()+ geom_density_2d(colour="black", size=2)

#Sort res by padj values
res <- res[order(res$padj), ]

#Merge tables with normalised counts
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

#load biomaRt to add other gene attributes in addition to Ensembl Gene ID
ensembl = biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

#Add gene symbols and biotypes
gene = getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"),values = resdata$Gene, mart = ensembl)
id <- base::match(resdata$Gene, gene$ensembl_gene_id)
resdata$Gene_name <- gene$external_gene_name[id]
resdata$Gene_biotype <- gene$gene_biotype[id]
resdata = resdata %>% relocate(Gene_name, .after = Gene)
resdata = resdata %>% relocate(Gene_biotype, .after = Gene_name)
head(resdata)

#output DESeq2 table
write.csv(resdata, file = "PCS_DESeq2.csv")

#Enhanced Volcano - |log2FC>1| & |padj<0.1| as cut off here

EnhancedVolcano(resdata,
                lab = as.character(resdata$Gene_name),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Recurrent vs Nonrecurret - No covariates',
                subtitle = NULL,
                selectLab = c('CD8B','PDCD1','GZMK','TOX','APOB'),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                xlim = c(-10,10),
                ylim = c(0,5),
                pCutoff = 0.1,
                FCcutoff = 1.0,
                pointSize = 1.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')
