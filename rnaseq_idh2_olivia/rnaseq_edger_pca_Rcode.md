# R Code for Differential Gene Expression and Principal Components Analysis

#### Preparing design

The design matrix contains value either 0 or 1 and will have one row for each sample and one column for each
combination of treatment conditions. Each sample (row) is assigned to one and one only
combination (column). 

```
R
library(edgeR)
source('~/svn_git/bioinformatics-cafe/trunk/makeTransparent.R')
source('~/svn_git/bioinformatics-cafe/trunk/intensityFilter.R')

targets<- read.table('data/exp_design.20140127.txt', header= TRUE, sep= '\t', row.names= 1)

## Rename columns to less verbose:
names(targets)[names(targets) == 'exp_design.hydroxygluterate_feeding']<- 'hdg'
names(targets)[names(targets) == 'exp_design.idh2_mutant']<- 'idh'
names(targets)[names(targets) == 'exp_design.transient_transfection']<- 'tt'
names(targets)[names(targets) == 'exp_design.photocaged']<- 'ph'
names(targets)[names(targets) == 'exp_design.light_time']<- 'lt' ## NB: Renamed "harvest_time" in django db
names(targets)[names(targets) == 'exp_design.cell_line']<- 'cl'
names(targets)[names(targets) == 'exp_design.amino_acid']<- 'aa'
names(targets)[names(targets) == 'exp_design.light']<- 'light'

targets$idh<- sub('_cell_line', '_cl', targets$idh)

## Remove cell_line column:
targets<- targets[, !names(targets) %in% 'cl']

## Defining each treatment combination as a group 
## See edgeRUserGuide 3.3.1
Group<- rep(NA, nrow(targets))
for(i in 1:nrow(targets)){
    g<- rep(NA, ncol(targets))
    for(j in 1:ncol(targets)){
        g[j]<- paste(names(targets)[j], targets[i,j], sep= '_')
    }
    Group[i]<- paste(g, collapse= '.')
}
Group<- as.factor(Group)

## For ease of reading: Table of contrats:
x<-  sort(as.character(unique(Group)))
xx<- cbind(x, as.matrix(do.call(rbind, strsplit(x, ".", fixed= TRUE))))
write.table(xx, file= 'design_conditions.txt', row.names= FALSE, col.names= FALSE, sep= '\t', quote= FALSE)

## Design matrix:
## --------------
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
row.names(design)<- row.names(targets)

## Check:
## ------
stopifnot(rowSums(design) == 1) ## Each sample has one and only one combination of conditions
colSums(design) ## Number of replicates per combination of conidtions: Must be 4 apart from hdg_no.idh_MUT.tt_yes.ph_no.light_no.lt_72h.aa_na which has 3 reps.
## ------

write.table(cbind(library_id= row.names(targets), design), 'design.txt', sep= '\t', col.names= TRUE, row.names= FALSE)
```

### Preparing data

Input file `htseq.matrix` is a tab delimited table with one row per gene and one
column per sample. Each entry is the raw count of reads assigned to each gene and
each sample. First column is gene name. This file might not be available in this
repository.

The file `hg19.gene_length.txt` has the length of each gene and it was produced
with [geneLengthFromGTF.py](https://github.com/dariober/bioinformatics-cafe/blob/master/geneLengthFromGTF.py)
using the appropriate reference gtf file. In this case we used `genes.gtf` from
[Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html) (Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf).

```
raw.data <- read.delim("data/htseq.matrix")
d <- raw.data[, 2:ncol(raw.data)] 
rownames(d) <- raw.data[, 1] 

## NOTE: Remove IDH2
# d<- d[which(row.names(d) != 'IDH2'),]

d <- DGEList(counts = d, group = Group)
dim(d)
cpm.d <- cpm(d) 

filter<- rowSums(cpm.d > 1) >=3
## filter[which(rownames(cpm.d)=='IDH2')]<- FALSE

d <- d[filter, ] 
d <- calcNormFactors(d)

# ----------------------------------------
# Gene expression
# ----------------------------------------
cpm.d <- cpm.d[filter, ] 
gene_length<- read.table('hg19.gene_length.txt', sep= '\t', col.names= c('gene_name', 'length'))
gene_length<- gene_length[gene_length$gene_name %in% row.names(cpm.d), ]
stopifnot(gene_length$gene_name ==  row.names(cpm.d))
rpkm.d<- rpkm(d, gene_length$length, log= TRUE)

write.table(cbind(row.names(rpkm.d), rpkm.d), 'gene_expression.rpkm.20140127.txt',col.names= TRUE, row.names= FALSE, sep= '\t', quote= FALSE)

## Grouping
pdf('hist_gene_expr.pdf')
hist(rpkm.d)
dev.off()
```

### Data exploration and principal components

```
## Dendrogram
hc<- hclust(dist(t(cpm.d)))

pdf('hclust.plot.pdf', width= 12/2.54, height= 26/2.54, pointsize= 8)
par(mar= c(2,1,4,6))
plot(as.dendrogram(hc), main= 'Hierarchical clustering\nof samples on genes counts per million', horiz= TRUE)
dev.off()

## PCA
pcaResult<-prcomp(t(rpkm.d[, 1: (which(colnames(rpkm.d) == 'osw112_65')-1) ]))
write.table(x= cbind(library_id= row.names(pcaResult$x), pcaResult$x), file= 'pcaResult.txt', quote= FALSE, sep= '\t',
    col.names= TRUE, row.names= FALSE)

pdf('pca.plot.cpm.pdf')
plot(pcaResult$x,
    main= 'Principal components of samples from photocage experiment',
    xlab= sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
    ylab= sprintf('PC2 (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
    type= 'n'
)
text(x= pcaResult$x[,1], y= pcaResult$x[,2], labels= rownames(pcaResult$x), cex= 0.5, col= c('red', 'blue', 'grey40'))
dev.off()

# PCA samples 1-48
# ----------------
pcaResult<-prcomp(t(rpkm.d[, 1: (which(colnames(rpkm.d) == 'osw095_48')) ]))
pdf('pca.plot.cpm_1-48.pdf')
plot(pcaResult$x,
    main= 'Principal components of samples from photocage experiment',
    xlab= sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
    ylab= sprintf('PC2 (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
    type= 'n'
)
text(x= pcaResult$x[,1], y= pcaResult$x[,2], labels= rownames(pcaResult$x), cex= 0.5, col= c('red', 'blue', 'grey40'))
dev.off()

## Exclude osw082 altogether
## -------------------------
dat<- rpkm.d[, 1: (which(colnames(rpkm.d) == 'osw095_48')) ]
dat<- dat[, which(colnames(dat) != 'osw082_35')]
pcaResult<-prcomp(t(dat))

pdf('pca.plot.cpm_1-48.wo_35.pdf')
plot(pcaResult$x,
    main= 'Principal components of samples from photocage experiment',
    xlab= sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
    ylab= sprintf('PC2 (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
    type= 'n'
)
text(x= pcaResult$x[,1], y= pcaResult$x[,2], labels= rownames(pcaResult$x), cex= 0.5, col= c('red', 'blue', 'grey40'))
dev.off()

## Zoom in to exclude osw082
## 
pdf('pca.plot.cpm_1-48.zoom.pdf')
plot(pcaResult$x,
    main= 'Principal components of samples from photocage experiment',
    xlab= sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
    ylab= sprintf('PC2 (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
    type= 'n',
    ylim= c(-40, 20)
)
text(x= pcaResult$x[,1], y= pcaResult$x[,2], labels= rownames(pcaResult$x), cex= 0.5, col= c('red', 'blue', 'grey40'))
dev.off()

# MDS plot
# --------
pdf('edger_mdsplot.pdf', width= 12/2.54, height= 12/2.54)
par(las= 1, cex= 0.8)
plotMDS(d, , xaxp= c(-2, 2, 8),  yaxp= c(-2, 2, 8), 
    col= c('blue', 'red'), 
    xlab= '', ylab= '',
    main= 'MDS plot of photocage experiment', 
    cex= 0.2)
grid()
dev.off()
```

### Gene Expression Analysis


```
# ---------------------------------------
# Estimating dispersion
# ---------------------------------------

d<- estimateGLMCommonDisp(d, design)
d<- estimateGLMTagwiseDisp(d, design)
fit<- glmFit(d, design)

# ---------------------------------------
# Contrasts
# ---------------------------------------
## MEMO: Contrasts are ordered alphanumerically
## E.g. `idh_yes - idh_no` is tested as  `idh_no - idh_yes` 
## hence +ve logFC means `idh_no > idh_yes`

## DrugvsPlacebo.2h = (Drug.2h-Drug.0h)-(Placebo.2h-Placebo.0h),

contrasts<- makeContrasts(
    x1_vs_2= hdg_no.idh_WT.tt_yes.ph_empty.light_no.lt_8h.aa_yes - hdg_no.idh_WT.tt_yes.ph_empty.light_yes.lt_8h.aa_yes,
    x7_vs_8= hdg_no.idh_WT.tt_yes.ph_yes.light_no.lt_8h.aa_yes - hdg_no.idh_MUT.tt_yes.ph_yes.light_yes.lt_8h.aa_yes,

    xA_vs_D= (hdg_no.idh_WT.tt_yes.ph_empty.light_no.lt_8h.aa_yes + hdg_no.idh_WT.tt_yes.ph_empty.light_yes.lt_8h.aa_yes) - 
             (hdg_no.idh_WT.tt_yes.ph_yes.light_no.lt_8h.aa_yes + hdg_no.idh_MUT.tt_yes.ph_yes.light_yes.lt_8h.aa_yes),

    x13_vs_14= hdg_no.idh_WT_WT.tt_yes.ph_no.light_no.lt_72h.aa_na - hdg_no.idh_MUT.tt_yes.ph_no.light_no.lt_72h.aa_na,
    levels= design
)

## Compute contrasts and write out DE tables
for (cntr in colnames(contrasts)[3]){
    lrt<- glmLRT(fit, contrast= contrasts[, cntr])
    detable<- topTags(lrt, n= nrow(d))$table
    detable$zScore<- localZ(detable$logCPM, detable$logFC, nbins= 20)

    pdf(paste('maplot', cntr, 'edgeR.pdf', sep= '.'))
    smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC', main= cntr)
    points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.01, 'grey', NA), pch= '.')
    points(x= detable$logCPM, y= detable$logFC, col= ifelse(abs(detable$zScore) > 2 & detable$FDR < 0.01, 'red', NA), pch= '.')
    legend('topright', pch= 19, col= c('grey', 'red'), legend= c(paste('N. FDR < 0.01=', sum(detable$FDR < 0.01)), paste('N. FDR < 0.01 & Z > 2=', sum(detable$FDR < 0.01 & abs(detable$zScore) > 2)) ))
    grid(col= 'grey50')
    dev.off()
    print(cntr)
    print(sum(detable$FDR < 0.01 & abs(detable$logFC) > log2(1.5)))

    write.table(cbind(gene_name= rownames(detable), detable), paste('expr_diff', cntr, 'edgeR.txt', sep= '.'), 
        row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
}
```