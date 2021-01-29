library(edgeR)
library(tidyverse)
library(openxlsx)
library(ggbiplot)
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/tschemic/Additional_Scripts/master/plot_cleanup.R", ssl.verifypeer = FALSE)))
#par(mfrow=c(1,1))

pval = 0.01 # set pvalue threshold (only relevant for plotting)
pval_adjust = "BH"  # set method for p-value adjustment for multiple testing (= FDR calculation)
cutoff = c(-1,1)  # set log2 fold change line in plots (e.g. "c(-1,1)" means lines will be drawn at 2-fold up and down)

targets <- readTargets(file = "Targets.txt") # this imports the data in the Targets.txt file = experiment setup

d <- readDGE(targets, header = FALSE) # this reads in the count data
# Filter low expression tags (cpm<1)
keep <- rowSums(cpm(d)> 1) >= 3 # this creates an R object with the genes/transcripts that have more than 1 counts per million - change the 1 if you want to modify the filter; the 3 at the end is the group size - change this if you have more or less replicates than 3
d <- d[keep,] # this filters out all genes/transcripts with less than 1 counts per million
d$samples$lib.size <- colSums(d$counts) # this calculates the library size again after the filtering step

# Normalization (TMM)
d <- calcNormFactors(d) # this calculates the normalization factors used during the diffrential expression analysis

# Estimating the dispersions and plot them
d <- estimateCommonDisp(d, verbose=FALSE) # this calculates the common dispersion over all genes and samples - is used during the diff. expr. anaylsis
d <- estimateTagwiseDisp(d) # this calculates the dispersion for each gene/transcript - is used during the diff. expr. anaylsis

pdf("disp.pdf")# change the name of the .pdf file; this opens a pdf file to store the next plot
plotBCV(d)
dev.off()# closes and saves the pdf file with the plot

# Differential expression, lines with # signs are just control commands to check the analysis

et140 <- exactTest(d, pair=c('CBS109','470140')) # this command performs the differential expression analysis
et147 <- exactTest(d, pair=c('CBS109','470147'))
et154 <- exactTest(d, pair=c('CBS109','470154'))

#detags <- rownames(topTags(et, n=20))
result140 <- summary(de140 <- decideTestsDGE(et140, p=pval, adjust=pval_adjust)) # this command gives the number of differentially expressed genes; add "lfc = 1" for logFC -1/1 as additional cutoff
result147 <- summary(de147 <- decideTestsDGE(et147, p=pval, adjust=pval_adjust))
result154 <- summary(de154 <- decideTestsDGE(et154, p=pval, adjust=pval_adjust))

detags140 <- rownames(d)[as.logical(de140)] # this command saves the row names of differentially expressed genes as "detags"
detags147 <- rownames(d)[as.logical(de147)]
detags154 <- rownames(d)[as.logical(de154)]

# generate smear plots (cloud of dots in log2CPM~log2FC axes)

pdf('470140.pdf')# change the name of the .pdf file; this opens a pdf file to store the next plot
plotSmear(et140, de.tags=detags140, xlab = "Average log2 counts per million", ylab = "log2 fold change") # this plots the analysis results; use "panel.first = NULL" to remove grid
abline(h = cutoff, col = "blue")
dev.off() # closes and saves the pdf file with the plot 

pdf('470147.pdf')# change the name of the .pdf file; this opens a pdf file to store the next plot
plotSmear(et147, de.tags=detags147, xlab = "Average log2 counts per million", ylab = "log2 fold change") # this plots the analysis results; use "panel.first = NULL" to remove grid
abline(h = cutoff, col = "blue")
dev.off() # closes and saves the pdf file with the plot 

pdf('470154.pdf')# change the name of the .pdf file; this opens a pdf file to store the next plot
plotSmear(et154, de.tags=detags154, xlab = "Average log2 counts per million", ylab = "log2 fold change") # this plots the analysis results; use "panel.first = NULL" to remove grid
abline(h = cutoff, col = "blue")
dev.off() # closes and saves the pdf file with the plot 

# Export results
res140 <- as.data.frame(topTags(et140, n=Inf)) # this saves the fold change table in "tt"
res147 <- as.data.frame(topTags(et147, n=Inf))
res154 <- as.data.frame(topTags(et154, n=Inf))

res <- merge(res140, res147, by=0)
res <- merge(res, res154, by.x=1, by.y=0)
colnames(res) <- c('geneID', 'logFC_470140', 'logCPM_470140', 'PValue_470140', 'FDR_470140',
 'logFC_470147', 'logCPM_470147', 'PValue_470147', 'FDR_470147',
 'logFC_470154', 'logCPM_470154', 'PValue_470154', 'FDR_470154')

write.xlsx(res, 'results.xlsx')

# generate volcano plots in (log2FC, -log10FDR) axes

pdf(file = "470140_FDR.pdf")
ggplot(res, aes(logFC_470140, -log10(FDR_470140))) + geom_point(alpha=0.4) + cleanup + 
geom_hline(yintercept = -log10(pval), colour='red') + geom_vline(xintercept = cutoff, colour='blue') + 
xlab('Log2 FC') + ylab('-Log10 FDR')
dev.off()

pdf(file = "470147_FDR.pdf")
ggplot(res, aes(logFC_470147, -log10(FDR_470147))) + geom_point(alpha=0.4) + cleanup + 
geom_hline(yintercept = -log10(pval), colour='red') + geom_vline(xintercept = cutoff, colour='blue') + 
xlab('Log2 FC') + ylab('-Log10 FDR')
dev.off()

pdf(file = "470154_FDR.pdf")
ggplot(res, aes(logFC_470154, -log10(FDR_470154))) + geom_point(alpha=0.4) + cleanup + 
geom_hline(yintercept = -log10(pval), colour='red') + geom_vline(xintercept = cutoff, colour='blue') + 
xlab('Log2 FC') + ylab('-Log10 FDR')
dev.off()


cpmill = as.data.frame(d$pseudo.counts) # this saves the counts per million table in a new R object called "pcounts_1h" - change this name as you want
colnames(cpmill) <- d$samples$group
write.xlsx(cpmill, 'cpmill.xlsx', rowNames=TRUE)

### PCA analysis

cpmill_transp <- t(cpmill)
cpmill_transp.pca <- prcomp(cpmill_transp, center = TRUE, scale. = TRUE)

plot <- ggbiplot(cpmill_transp.pca, var.axes = FALSE)
plot + cleanup + geom_text(aes(label=d$samples$description))
pdf(file = "PCA.pdf")
plot + cleanup + geom_text(aes(label=d$samples$description))
dev.off()



### Pairwise differential expression analysis of 470140 vs. 470154
library(edgeR)
library(tidyverse)
library(readxl)
library(RCurl)
library(ggsci)
library(ggbiplot)
eval(parse(text = getURL("https://raw.githubusercontent.com/tschemic/Additional_Scripts/master/plot_cleanup.R", ssl.verifypeer = FALSE)))


pval = 0.01 # set pvalue threshold (only relevant for plotting)
pval_adjust = "BH"  # set method for p-value adjustment for multiple testing (= FDR calculation)
cutoff = c(-1,1)  # set log2 fold change line in plots (e.g. "c(-1,1)" means lines will be drawn at 2-fold up and down)

targets <- readTargets(file = "Targets.txt") # this imports the data in the Targets.txt file = experiment setup
targets <- targets[c(1:3,7:9),]

d <- readDGE(targets, header = FALSE) # this reads in the count data
# Filter low expression tags (cpm<1)
keep <- rowSums(cpm(d)> 1) >= 3 # this creates an R object with the genes/transcripts that have more than 1 counts per million - change the 1 if you want to modify the filter; the 3 at the end is the group size - change this if you have more or less replicates than 3
d <- d[keep,] # this filters out all genes/transcripts with less than 1 counts per million
d$samples$lib.size <- colSums(d$counts) # this calculates the library size again after the filtering step

# Normalization (TMM)
d <- calcNormFactors(d) # this calculates the normalization factors used during the diffrential expression analysis

# Estimating the dispersions and plot them
d <- estimateCommonDisp(d, verbose=FALSE) # this calculates the common dispersion over all genes and samples - is used during the diff. expr. anaylsis
d <- estimateTagwiseDisp(d) # this calculates the dispersion for each gene/transcript - is used during the diff. expr. anaylsis

pdf("disp_140vs154.pdf")# change the name of the .pdf file; this opens a pdf file to store the next plot
plotBCV(d)
dev.off()# closes and saves the pdf file with the plot

# Differential expression, lines with # signs are just control commands to check the analysis

et <- exactTest(d, pair=c('470140','470154')) # this command performs the differential expression analysis

#detags <- rownames(topTags(et, n=20))
result <- summary(de <- decideTestsDGE(et, p=pval, adjust=pval_adjust)) # this command gives the number of differentially expressed genes; add "lfc = 1" for logFC -1/1 as additional cutoff

detags <- rownames(d)[as.logical(de)] # this command saves the row names of differentially expressed genes as "detags"

# generate smear plots (cloud of dots in log2CPM~log2FC axes)

pdf('470154vs470140.pdf')# change the name of the .pdf file; this opens a pdf file to store the next plot
plotSmear(et, de.tags=detags, xlab = "Average log2 counts per million", ylab = "log2 fold change") # this plots the analysis results; use "panel.first = NULL" to remove grid
abline(h = cutoff, col = "blue")
dev.off() # closes and saves the pdf file with the plot 

# Export results
res <- as.data.frame(topTags(et, n=Inf)) # this saves the fold change table in "tt"
rownames(res) <- gsub(pattern = "gene-", replacement = "", x = rownames(res))

write.table(res, "154vs140_results.tsv", sep = "\t")

# Add C. albicans orthologs

CaOrth <- read_excel("../required_files/Ortho.xlsx")
CaOrth <- CaOrth[,c(1,2,6)]
CaOrth <- unique(CaOrth)

resCa <- merge(res, CaOrth, by.x = 0, by.y = 2)

IDtable <- read_tsv("../required_files/A22_A19_names.txt")
IDtable$GeneName <- ifelse(is.na(IDtable$GeneName), IDtable$CGDID, IDtable$GeneName)

resNames <- merge(resCa, IDtable, by.x = 6, by.y = 1, all.x = TRUE)
write.table(resNames[,c(1,2,10,8,3,6)], "154vs140_results_wNames.tsv", sep = "\t", row.names = FALSE)

# generate volcano plots in (log2FC, -log10FDR) axes

# Correct FDR = 0 values to not get Inf after log transformation
resCa$FDR[resCa$FDR == 0] <- 1e-200

pdf(file = "470154vs470140_VulcPlot.pdf", height = 5)
ggplot(resCa, aes(logFC, -log10(FDR))) + geom_point(alpha=0.4) + cleanup + 
  geom_hline(yintercept = -log10(pval), colour='red') + geom_vline(xintercept = cutoff, colour='blue') + 
  xlab('Fold Change (log2)') + ylab('Adjusted P-Value (-log10)') +
  geom_text(aes(label=ifelse((resCa$logFC > 1.5 | resCa$logFC < -1.5 | -log10(resCa$FDR) > 40), 
                             resCa$ANNOT_GENE_NAME, ""), vjust=1.5), size = 3)

dev.off()


cpmill = as.data.frame(d$pseudo.counts) # this saves the counts per million table in a new R object called "pcounts_1h" - change this name as you want
colnames(cpmill) <- d$samples$description
rownames(cpmill) <- gsub(pattern = "gene-", replacement = "", x = rownames(cpmill))

write.table(cpmill, '154vs140_cpmill.tsv', sep = "\t")

### PCA analysis

cpmill_transp <- data.frame(t(cpmill))
colnames(cpmill_transp) <- rownames(cpmill)
cpmill_transp$group <- gsub(pattern = "_rep[1-3]", replacement = "", x = rownames(cpmill_transp))
cpmill_transp.pca <- prcomp(dplyr::select(cpmill_transp, -group), scale. = TRUE)

plot <- ggbiplot(cpmill_transp.pca, var.axes = FALSE, groups = cpmill_transp$group)
#plot + cleanup + geom_text(aes(label=d$samples$description))
plot2 <- plot + cleanup + theme(legend.key = element_blank()) +
  xlab("PC1 (52.4% explained variance)") + 
  ylab("PC2 (26.8% explained variance)") + 
  scale_color_brewer(palette = "Dark2", name = expression(paste(italic("C. auris"), " strains")), labels = c("Cau470140", "Cau470154"))
plot2

pdf(file = "470154vs470140_PCA.pdf", height = 5)
plot2
dev.off()

png(file = "470154vs470140_PCA.png", res = 300, width = 1920, height = 1440)
plot2
dev.off()


