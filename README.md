---
title: "Project Stat597A. n\ Bioinformatics for High Throughput Experiments"
author: "Hannier Pulido & Yong Jung"
date: "April 22, 2014"
output: pdf_document
---

# **Objective**
The aim of this project is to identify a gene signature for the four different types of medulloblastoma described by (Robinson et al., 2012). The analysis reveals a list of genes whose expression is correlated with the medulloblastoma subgroups

# **Introduction**
In order to classify the samples in the four different types of medulloblastoma, we used the expression profiling by array available at the Gene Expression Omnibus website with the accession number GSE37418. A total of 54,675 genes from 76 samples were analyzed using Bioconductor and Random Forest (Breiman, 2001) for classification. 

(Develop more introduction with the results obtained from the annotation and gene ontology… what did we discover? What is the gene signature for each of the four groups?)

# **Step 1. Load libraries and CEL files**
```{r Step 1 Install packages, warning=FALSE, message=FALSE}
source("http://www.bioconductor.org/biocLite.R")
biocLite("simpleaffy")
biocLite("limma")
biocLite("affy")
biocLite("affy")
biocLite("hexbin")
biocLite("annotate")
biocLite("hgu133plus2.db") 
biocLite("GOstats")
biocLite("GO.db")
biocLite("genefilter")

# Install microarray analysis packages
require("simpleaffy")
require("limma")
require("affy")
require("hexbin")
require("annotate")
require("hgu133plus2.db") # Affymetrix Human Genome U133 Plus 2.0 Array annotation datar
require("genefilter")
require("simpleaffy")
require("GOstats")
require("GO.db")
require("genefilter")

# Install Random Forest packages
library(RRF)
library(randomForest)
library(varSelRF)
library(lattice)
library(Boruta)

# Other packages
library(scatterplot3d)
library(mclust)
library(gplots)
library(graphics)
library(ClassDiscovery)
library(lattice)
library(xtable)

# Upload the metadata describing samples in the experiment with information of subroup, age, sex, ethnic and Mstage 
phenoData <- read.AnnotatedDataFrame("annotation.txt", header = TRUE, sep = "\t", row.names = 1)

# Load raw CEL files using ReadAffy ("affy")
CELdata <- ReadAffy(celfile.path = 'MicroarrayData', phenoData = phenoData)
# cancer.rma <- rma(CELdata) # Optional step for background correction using RMA
```
https://github.com/hannierpulido/Project.wiki.git
# **Step 2. Quality assessment of raw data**
We assess the quality of the raw data before correcting and normalizing the data for further analysis. This analysis allow us to detect abnormal arrays on the data

## **2.1 Boxplots of raw data**
```{r Step 2 Boxplots of raw data, warning=FALSE, message=FALSE}
# Generates a vector for colours
pdata <- pData(CELdata)

color <- NULL
for(i in 1:length(pdata$subgroup)) {
  clas <- pdata$subgroup[i]
  if(clas == "G3") {
    color = c(color, "green3")
  } else if(clas == "G4") {
    color = c(color, "cyan2")
  } else if(clas == "SHH") {
    color = c(color, "orange")
  } else if(clas == "SHH_OUTLIER") {
    color = c(color, "yellow")
  } else if(clas == "U") {
    color = c(color, "magenta")
  } else if(clas == "WNT") {
    color = c(color, "firebrick1")
  } else {
    color = c(color, "white")
  }
}

# Creates a vector with the name of genes extracted from the CELdata object
D <- data.frame(sampleNames(CELdata))
D$genLabel <- lapply(strsplit(as.character(D$sampleNames.CELdata.), "\\_"), "[", 1) 

par(mar = c(5.1, 4.1, 4.1, 2.1))
boxplot(CELdata, col = color, las = 3, xaxt = "n", main = "Boxplot of raw data")
axis(1, at = 1:76, labels = D$genLabel, las = 3, cex.axis = 0.6)
legend(-2, 16, horiz = TRUE, c("G3", "G4", "SHH", "SHH_OUTLIER", "U", "WNT"), 
       fill = c("green3", "cyan2", "orange", "yellow","magenta", "firebrick1") )
```

## **2.2 Density plots of raw data**
```{r Step 2 Density plots of raw data, warning=FALSE, message=FALSE}
hist(CELdata, type = "l", col = color, main = "Histogram of raw data")
legend("topright", c("G3", "G4", "SHH", "SHH_OUTLIER", "U", "WNT"), 
       fill = c("green3", "cyan2", "orange", "yellow","magenta", "firebrick1") )
```

## **2.3 RNA degradation Plot**
In this plot, we are looking if any arrays are really different from the others. 
```{r step 2 RNA degradation plot, warning=FALSE, message=FALSE}
RNA.deg <- AffyRNAdeg(CELdata)
plotAffyRNAdeg(RNA.deg, col = color)
legend("topleft", c("G3", "G4", "SHH", "SHH_OUTLIER", "U", "WNT"), 
       fill = c("green3", "cyan2", "orange", "yellow","magenta", "firebrick1") )
```
```
Conclusion: data needs to be normalized before further analysis

# **3. Preprocess: Background correction and normalization**
The expression data were filtered to remove any probe sets that failed to show significat variation in expression across the data set (Absent genes.
Then, we used the function **rma** (*affy*) to perform RMA background correction and quantile normalization.   
Since we are only interested in the four medulloblastoma groups originally described by the authors, the unidentified groups and one outlier (labeled by the author) were not considered for further analysis.  

These two filtered steps reduced the complexity of the data from 76 samples and 54,675 probe sets to 73 samples and  47,330 probe sets
```{r Step 3 Preprocess, warning=FALSE, message=FALSE}
# Get the ExpressionSet object using RMA background correction and quantile normalization
cancer.rma <- rma(CELdata)
dim(cancer.rma)

#===============================================================================
#        Remove absent genes (No variation across any probe sets)
#===============================================================================
# The expression profile data is improved by filtering out the absent genes identified by *mas5calls*
# Mas5calls returns the presence/absence states of all genes in each array

# Performs the Wilcoxon signed rank-based gene expression presence/absence
eset.mas5 <- mas5calls(CELdata)

# Creates a function to remove the absent genes
mascallsfilter <- function(cutoff = "A", number){
          function(x){
          sum((x) == (cutoff)) != number
          }
    } 

# Filtering the absent genes       
m5cfun <- mascallsfilter(number = dim(eset.mas5)[2])     
mcfun <- filterfun(m5cfun)
mcsub <- genefilter(eset.mas5, mcfun)
sum(mcsub) # 47330
 
# subset of the cancer.rma ExpressionSet (it no longer has absent genes)
cancer.sub1 <- cancer.rma[mcsub,]
dim(cancer.sub1)        # 47330    76
# After the gene filtering we reduced the dimensions of the database from 54,675 to 47,330 genes

#===============================================================================
#        Remove "SHH_OUTLIER" & "U"
#===============================================================================
cancer.sub <- get.array.subset.affybatch(cancer.sub1, "subgroup", c("G3", "G4", "SHH", "WNT"))
dim(cancer.sub)

# Creates a vector with the name of genes extracted from the cancer.sub object
C <- data.frame(sampleNames(cancer.sub))
C$genLabel <- lapply(strsplit(as.character(C$sampleNames.cancer.sub.), "\\_"), "[", 1) 

# Change the sample names in the cancer.sub object
sampleNames(cancer.sub) <- as.character(C$genLabel)

# We will be working with **cancer.sub** from now on
```

# **Step 4. Preliminary analysis**
A clustering analysis of k-means using 2 and 4 groups were performed to explore the association between the resulting cluster sets and the known cancer subgroups.

We used Fisher's exact test to analyze the difference between the 2 and 4 clustering and decided to use 4 groups
```{r Step 4 Kmeans clustering, warning=FALSE, message=FALSE}
# Extracts metadata from cancer.sub expressionSet
pdata <- pData(cancer.sub)
expr <- data.frame(exprs(cancer.sub)) # Expression data
label <- as.vector(pData(cancer.sub)$subgroup) # character vector for the subgroups

# kmeans for 2 clusters
cluster_2 <- kmeans(t(expr), iter.max = 1000, 2)$cluster

# Combine metadata with cluster for 2 kmeans
res_2 <- cbind(pdata, cluster_2)
# write.table(res_2, file = "Cluster_2_Annotation.txt", sep = '\t')
table(res_2$subgroup, res_2$cluster)

# kmeans for 4 clusters
cluster_4 <- kmeans(t(expr), iter.max = 1000, 4)$cluster

# Combine metadata with cluster for 4 kmeans
res_4 <- cbind(pdata, cluster_4)
# write.table(res_4, file = "Cluster_4_Annotation.txt", sep = '\t')
table(res_4$subgroup, res_4$cluster)

# Fisher's exact test
fisher.test(table(res_2$subgroup, res_2$cluster_2), workspace=2e9)
fisher.test(table(res_4$subgroup, res_4$cluster_4), workspace=2e9)

#===============================================================================
#        MDS plot for k-means
#===============================================================================

# Generates a new vector for colors
color <- NULL
for(i in 1:length(pdata$subgroup)) {
  clas <- pdata$subgroup[i]
  if(clas == "G3") {
    color = c(color, "green3")
  } else if(clas == "G4") {
    color = c(color, "cyan2")
  } else if(clas == "SHH") {
    color = c(color, "orange")
  } else if(clas == "WNT") {
    color = c(color, "firebrick1")
  } else {
    color = c(color, "white")
  }
}

# Generates the distance matrix
dist_res <- dist(t(expr))
mds<-cmdscale(dist_res) # Principal coordinates analysis for the distance matrix
x<-mds[,1]
y<-mds[,2]

plot(x, y, type = 'n', main = "MDS plot according to a k-means result")
text(x, y, label, col = color, cex = 1, font = 2)
```

# **Step 5. Quality assessment of filtered data**

## **5.1 Boxplots of filtered data**
```{r Step 5 Boxplot filtered data, warning=FALSE, message=FALSE}
boxplot(exprs(cancer.sub), col = color, las = 3, xaxt = "n", main = "Boxplot of filtered data")
axis(1, at = 1:73, labels = C$genLabel, las = 3, cex.axis = 0.6)
legend(20, 15, horiz = TRUE, bg = "white", c("G3", "G4", "SHH", "WNT"), 
       fill = c("green3", "cyan2", "orange", "firebrick1") )
```

## **5.2 Density plots of filtered data**
```{r Step 5 Density plots, , warning=FALSE, message=FALSE}
plot(density(expr[, 1]), col = color[[1]], ylim = c(0, 0.25), xlab = "", 
     main = "Histogram of filtered data")
for(i in 2:ncol(expr)) lines(density(expr[, i]), col = color[i])
legend("topright", c("G3", "G4", "SHH", "WNT"), 
       fill = c("green3", "cyan2", "orange", "firebrick1") )
```

## **5.3 MVA plots for filtered data**
```{r step 5,  warning=FALSE, message=FALSE}
# Obtain the number of samples for the different groups
which(pdata$subgroup == "G3") # Returns the number of columns with samples G3
which(pdata$subgroup == "G4") # Returns the number of columns with samples G4
which(pdata$subgroup == "SHH") # Returns the number of columns with samples SHH
which(pdata$subgroup == "WNT") # Returns the number of columns with samples WNT

# MVA plot for four arrays for the group G3
mva.pairs(expr[, c(9, 23, 46, 61)], col = color, cex = 1, 
          main = "MVA plot for 4 arrays of G3") 

# MVA plot for four arrays for the group G4
mva.pairs(exprs(cancer.sub)[, c(1, 33, 27, 56)], col = color, cex = 1, 
          main = "MVA plot for 4 arrays of G4") 

# MVA plot for four arrays for the group SHH
mva.pairs(exprs(cancer.sub)[, c(5, 31, 42, 71)], col = color, cex = 1, 
          main = "MVA plot for 4 arrays of SHH") 

# MVA plot for four arrays for the group WNT
mva.pairs(exprs(cancer.sub)[, c(3, 47, 62, 64)], col = color, cex = 1,
          main = "MVA plot for 4 arrays of WNT") 
```

# **Step 6. Exploratory analysis of the pre-processed data**
We used an unsupervised 2-D hierarchical clustering (HCA) to assess the major grouping of the four types of medulloblastoma based only on their gene expression profiles.
This exploratory analysis supports the K-means clustering previously performed in step

## **6.1 Hierarchical clustering**
```{r step 6 HCA,  warning=FALSE, message=FALSE}
A <- t(expr)
dim(A)

label <- as.vector(pData(cancer.sub)$subgroup)

# Computes the distance matrix
my.dist <- function(x) dist(x, method = "euclidean")

# Hierarchical clustering on the filtered data
my.hclust <- function(d) hclust(d, method = "ward")

# Generates Hierarchical cluster analysis dendrogram
d <- as.dist(dist(A, method = "euclidean"))
fit.h <- hclust(d, method = "ward")
plotColoredClusters(fit.h, labs = label, cols = color, cex = 0.8, 
                    xlab = "", sub = "", main = "Hierarchical clustering")
```

## **6.2 PCA**
```{r step 6 PCA,  warning=FALSE, message=FALSE}
pca <- prcomp(A, scale = T)
summary(pca)
scatterplot3d(pca$x[, 1:3], pch = 20, cex.symbols = 2, color = color, 
              main = "PCA for the filtered genes")
legend(5.5, 3, c("G3", "G4", "SHH", "WNT"), 
       fill = c("green3", "cyan2", "orange", "firebrick1") )
```
Both, the PCA and HCA imply that the samples are clustered in four different groups, which agrees with the original subgrouping


# **Step 7. Random Forest**
The genes selected as important features by the random forest are the genes that support the prediction of the outcome, i.e. the separation of the four groups.
We can subset the data to only inlclude the most important variables and use that with another model (like ANOVA, GLM).  
The measure of importance means that if the variable is not important (the null hypothesis), then rearranging the values of that variable will not degrade prediction accuracy.
```{r Step 7 RF,  warning=FALSE, message=FALSE}
subgroup <- factor(pdata$subgroup)
# Tune Random Forest to get the right mtry for the lowest OOB
set.seed(666)
tuneRRF(t(exprs(cancer.sub)), factor(subgroup)) # for mtry = 217, the OOB is 9.59%

# Computes the random forest
set.seed(666)
rf <- randomForest(t(exprs(cancer.sub)), factor(subgroup), mtry = 217, ntree = 500,
                   importance = TRUE, proximity = TRUE) 

# Plot the error rage (y-axis) vs the number of trees
colRF <- c("black", "green3", "cyan2", "orange", "firebrick1")
plot(rf, col = colRF, main = "random Forest error rate")
legend("topright", colnames(rf$err.rate), col = colRF, fill = colRF)

# Plot the clustering resulting from Random Forest
randomForest:::MDSplot(rf, factor(subgroup), main = "Clustering plot for Random Forest",
                       palette = c("green3", "cyan2", "orange", "firebrick1"), cex = 2)
legend("topleft", horiz = TRUE, bg = "white", c("G3", "G4", "SHH", "WNT"), 
       fill = c("green3", "cyan2", "orange", "firebrick1") )
```

# Step 8. Gene Signature detection
After fitting the Random Forest model, two approaches were used for the variable selection:
 
- **VarSelRF** Uses the OOB error as minimization criterion. It carry out variable elimination of the least important variables
 
- **Boruta** inds all relevant genes on the learning algorithm (RF) that are “biologically” related to the medulloblastoma groups. Uses a p-value cut-off of 0.05
```{r step 8 Gene Signature detection, warning=FALSE, message=FALSE}
# Variable selection from random forest using OOB error
# This variable selection is performed using the OOB error as minimization criterion. 
set.seed(666)
rf.vs1 <- varSelRF(t(exprs(cancer.sub)), factor(subgroup))

# Get the most important genes according to VarSelRF
Best.vars <- data.frame(rf.vs1$selected.vars)
colnames(Best.vars) <- "Best"
head(Best.vars)
dim(Best.vars) # 58 genes are selected by VarSelRF

#===============================================================================
#        Plot of the importance measure for the first 500 genes
#===============================================================================

# Getting the variable importance values (MeanDecreaseAccuracy from varSelRF)
imp <- rf.vs1$firstForest$importance
imp2 <- sort(imp[, 5], decreasing = TRUE) 

var.imp <- dotplot(sort(imp2[1:500]), 
                   main = "Important genes for classifying four medulloblastoma groups",
                   scales = list(y = list(cex = 0.5)),
                   xlab = list(label = "Mean Decrease of Accuracy", cex = 1))

update(var.imp, panel = function(...) {
  panel.abline(h = 380, v = 0.0005, lty = "longdash", col = "red", lwd = 2)
  panel.xyplot(...)
})

# This plot is showing the first 500 most important genes ranked according to the 
# Mean Decrease of Accuracy (measure of importance).
# The red lines indicate the elbow of the curve. Genes above this threshold can 
# be considered as the most important predictors for the RF. 
# However, to avoid subjective thresholds, we used VarSelRF and boruta to select
# the most important gens.

#===============================================================================
#        variable selection using boruta
#===============================================================================
# Variable selection using Boruta package (Caution: Execution time is 1:45 hours in a 16GB PC)
set.seed(666)
Bor <- Boruta(t(expr), factor(subgroup), pValue = 0.05)

# Extract names of the selected attributes
Best.Bor <- data.frame(getSelectedAttributes(Bor, withTentative = FALSE))
colnames(Best.Bor) <- "Best"

# For some reason, an "X" is prefixed before the name of the gene.
# This code should eliminate the "X"
Best.Bor$Best <- lapply(strsplit(as.character(Best.Bor$Best), "\\X"), "[", 2) 
dim(Best.Bor) # 66 genes were selected by boruta
head(Best.Bor)

# Venn Diagram showing agreement between both variable selection techniques
venn(list(VarSelRF = Best.vars, Boruta = Best.Bor)) 
# 58 genes were selected by VarSelRF
# 66 genes were selected by Boruta
# 39 genes were common for both approaches (interSet)
# 85 genes is the union between the results from both approaches (unionSet)
```

## **8.1 Export the results of Random Forest**
```{r Step 7 Export results from RF,  warning=FALSE, message=FALSE}
# Merge the genes from varselRF and boruta to get the genes in common (interSet)
Best.Bor$Best <- as.character(Best.Bor$Best)
Best.vars$Best <- as.character(Best.vars$Best)
interSet<- merge(Best.Bor, Best.vars)
dim(interSet)
head(interSet)

# Export the bestGenes containing 39 genes in common
write.csv(interSet, file = "interSet.csv", row.names = FALSE)

# Combine the best genes from Boruta and varSelRF and drop the duplicates (unionSet)
unionSet <- c(Best.Bor$Best, Best.vars$Best)
length(unionSet)
unionSet <- data.frame(unique(unionSet))
colnames(unionSet) <- "best"
dim(unionSet)
head(unionSet)

# Export the bestGenes containing 85 genes among the two approaches (unionSet)
write.csv(unionSet, file = "unionSet.csv", row.names = FALSE)

Best81 <- read.csv("bestGenesPrev.csv")

venn(list(Previous = Best81, current = bestGenes))
```


# **Step 9 Heatmap of the best genes**
1. Extract expression set
2. Match names of genes in expression set with names in the bestGenes dataframe
3. Do the heatmap
```{r step 9 Heatmap of best genes, warning=FALSE, message=FALSE}
# Generates heatmap of the most important genes

bestGenes <- read.csv("bestGenes.csv")
B <- exprs(cancer.sub)
dim(B)

keep <- subset(B, rownames(B) %in% bestGenes[,1])
dim(keep)

bt <- t(keep)

hr <- heatmap.2(keep, labCol = label, scale = "row", col=greenred(75), ColSideColors = color, key = TRUE, symkey = FALSE, density.info = "none", trace = "none", main = "Heatmap for most important genes")

hr <- heatmap.2(keep, labCol = label, scale = "row", col=greenred(75), key = TRUE, symkey = FALSE, density.info = "none", trace = "none", main = "Heatmap for most important genes")

hr$rowInd
hr$colInd


# Extract the genes ID according to the two gene sets in the cluster of the heatmap
keep.names <- rownames(keep)
ID <- hr$rowInd
geneSets <- keep.names[ID]
geneSetsID <- data.frame(as.numeric(ID), geneSets)
colnames(geneSetsID) <- c("ID", "gene")
write.csv(geneSetsID, file = "geneSetsID.csv", row.names = FALSE)


####

####
# Heatmap of the best common genes

bestGenes.common <- read.csv("bestGenesCommon.csv")
keep2 <- subset(B, rownames(B) %in% bestGenes.common[,1])
dim(keep2)

hr2 <- heatmap.2(keep2, labCol = label, scale = "row", col=greenred(75), ColSideColors = color, key = TRUE, symkey = FALSE, density.info = "none", trace = "none", main = "Heatmap for most important genes \n common genes between VarSelRF and Boruta")

hr2$rowInd
hr2$colInd

# Extract the genes ID according to the two gene sets in the cluster of the heatmap
keep2.names <- rownames(keep2)
ID2 <- hr2$rowInd
geneSetsCommon <- keep.names[ID2]
geneSetsIDCommon <- data.frame(as.numeric(ID2), geneSetsCommon)
colnames(geneSetsIDCommon) <- c("ID", "gene")
write.csv(geneSetsIDCommon, file = "geneSetsIDCommon.csv", row.names = FALSE)


hc.rows<- hclust(dist(t(bt)))
)

ct<- cutree(hc.rows, h=10) # it gives me 6 groups

rect.hclust(hc.rows, h=10) # draw red rectangles to mark the subgroups

table(ct)

tableclust<- data.frame(t(bt),ct)


hc <- as.hclust( hm$rowDendrogram )
cutree( hc, h=10 )

########################################################
########################################################
######## Venn Diagram for initial 81 genes selected (This is to show the overlap of boruta and varselrf from the first run, I didn't set seeds at that time, so the results are different from set.seed(666). DON'T INCLUDE THIS SECTION IN THE FINAL REPORT####
########################################################
########################################################

par(mar=c(5.1,4.1,4.1,2.1)
par(mar=c(2,2,2,2))
par(oma=c(2,2,2,2))
frame()
draw.pairwise.venn(56, 34, 29, category = c("Boruta", "VarSelRF"), scaled = FALSE, cex = 3, cat.cex = 3, cat.pos = 0, fill = c("cyan", "orange"))

```


# Step 10. Gene Ontology
```{r step 10 Gene Ontology}
#GO enrichment

ls("package:hgu133plus2.db")
head(as.list(hgu133plus2GO))


interSet <- c("1559149_at", "1561341_at", "1562828_at", "201209_at", "201416_at", "202719_s_at", "203576_at", "203607_at", "205637_s_at", "205667_at", "212288_at", "212774_at", "215047_at", "217626_at", "221798_x_at", "225512_at", "225955_at", "226049_at", "227785_at", "231325_at", "231631_at", "231776_at", "233546_at", "235004_at", "238850_at", "241881_at", "243428_at", "243435_at", "244370_at")

unionSet <- c("1555392_at", "1559149_at", "1561341_at", "1562828_at", "200715_x_at", "200903_s_at", "201209_at", "201416_at", "202719_s_at", "203576_at", "203607_at", "204457_s_at", "204546_at", "205637_s_at", "205667_at", "206017_at", "206456_at", "212288_at", "212774_at", "215028_at", "215047_at", "215531_s_at", "217626_at", "218479_s_at", "218736_s_at", "220111_s_at", "221798_x_at", "222725_s_at", "224539_s_at", "224955_at", "225381_at", "225512_at", "225955_at", "226049_at", "226796_at", "227439_at", "227785_at", "228665_at", "230902_at", "231325_at", "231631_at", "231776_at", "232003_at", "232069_at", "232113_at", "232286_at", "233546_at", "235004_at", "235044_at", "238850_at", "241411_at", "241881_at", "242527_at", "243428_at", "243435_at", "244370_at", "200807_s_at", "200845_s_at", "202260_s_at", "204105_s_at", "205378_s_at", "205493_s_at", "206186_at", "208692_at", "210271_at", "212713_at", "214925_s_at", "215290_at", "216132_at", "218394_at", "221883_at", "221923_s_at", "224092_at", "227290_at", "227441_s_at", "230019_s_at", "230765_at", "241505_at", "242959_at", "243521_at", "244076_at")


a <- interSet

probeIds <- rownames(rat.rma)
GOannot = mget(probeIds, hgu133plus2GO)
notEmpty = function(x) {
  if (length(x) == 1 && is.na(x))
  	FALSE
	else TRUE
}

haveGo = sapply(GOannot,notEmpty)
GoFeatures = probeIds[haveGo]

entrezIds <- mget(GoFeatures, hgu133plus2ENTREZID)
haveEntrezId <- sapply(entrezIds, function(x) !is.na(x))
entrezFeatures = GoFeatures[haveEntrezId]

require(genefilter)
EntrezSelect = names(findLargest(a, rnorm(length(a)), "hgu133plus2.db"))
EntrezUniverse = names(findLargest(entrezFeatures, rnorm(length(entrezFeatures)), "hgu133plus2.db"))


InteractBPover = new("GOHyperGParams", geneIds = EntrezSelect, universeGeneIds = EntrezUniverse, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
InteractSigBPover = hyperGTest(InteractBPover)

summaryTab = summary(InteractSigBPover)
summaryTab[summaryTab$P <= 0.05, ]

#get("GO:0006550", GOTERM)

interRes <- summaryTab # inter_GO.txt



a <- unionSet

probeIds <- rownames(rat.rma)
GOannot = mget(probeIds, hgu133plus2GO)
notEmpty = function(x) {
	if (length(x) == 1 && is.na(x))
		FALSE
	else TRUE
}

haveGo = sapply(GOannot,notEmpty)
GoFeatures = probeIds[haveGo]

entrezIds <- mget(GoFeatures, hgu133plus2ENTREZID)
haveEntrezId <- sapply(entrezIds, function(x) !is.na(x))
entrezFeatures = GoFeatures[haveEntrezId]

require(genefilter)
EntrezSelect = names(findLargest(a, rnorm(length(a)), "hgu133plus2.db"))
EntrezUniverse = names(findLargest(entrezFeatures, rnorm(length(entrezFeatures)), "hgu133plus2.db"))


InteractBPover = new("GOHyperGParams", geneIds = EntrezSelect, universeGeneIds = EntrezUniverse, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
InteractSigBPover = hyperGTest(InteractBPover)

summaryTab = summary(InteractSigBPover)
summaryTab[summaryTab$P <= 0.05, ]

#get("GO:0006550", GOTERM)

unionRes <- summaryTab # union_GO.txt
write.table(unionRes, file = "union_GO.txt", sep = '\t')

#After heatmap with unionSet

unionSet_1 <- c("206456_at", "218736_s_at", "220111_s_at", "226796_at", "241505_at", "218479_s_at", "224539_s_at", "242527_at", "215047_at", "215531_s_at", "241881_at", "222725_s_at", "232069_at", "235044_at", "228665_at", "212713_at", "243428_at", "203576_at", "205667_at", "243435_at", "215028_at", "204457_s_at", "232286_at", "241411_at", "201209_at", "200903_s_at", "221923_s_at", "230765_at", "1559149_at", "235004_at", "206186_at", "242959_at", "227441_s_at", "227290_at")


unionSet_2 <- c("232003_at", "204546_at", "231631_at", "217626_at", "1555392_at", "206017_at", "216132_at", "205378_s_at", "202719_s_at", "244076_at", "243521_at", "225955_at", "233546_at", "204105_s_at", "200845_s_at", "224955_at", "238850_at", "205493_s_at", "212288_at", "202260_s_at", "212774_at", "201416_at", "200715_x_at", "200807_s_at", "208692_at", "221798_x_at", "231776_at", "225381_at", "232113_at", "230902_at", "215290_at", "1561341_at", "1562828_at", "231325_at", "224092_at", "244370_at", "205637_s_at", "218394_at", "210271_at", "230019_s_at", "227785_at", "214925_s_at", "227439_at", "203607_at", "225512_at", "226049_at", "221883_at")

a <- unionSet_1

probeIds <- rownames(rat.rma)
GOannot = mget(probeIds, hgu133plus2GO)
notEmpty = function(x) {
  if (length(x) == 1 && is.na(x))
		FALSE
	else TRUE
}

haveGo = sapply(GOannot,notEmpty)
GoFeatures = probeIds[haveGo]

entrezIds <- mget(GoFeatures, hgu133plus2ENTREZID)
haveEntrezId <- sapply(entrezIds, function(x) !is.na(x))
entrezFeatures = GoFeatures[haveEntrezId]

require(genefilter)
EntrezSelect = names(findLargest(a, rnorm(length(a)), "hgu133plus2.db"))
EntrezUniverse = names(findLargest(entrezFeatures, rnorm(length(entrezFeatures)), "hgu133plus2.db"))


InteractBPover = new("GOHyperGParams", geneIds = EntrezSelect, universeGeneIds = EntrezUniverse, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
InteractSigBPover = hyperGTest(InteractBPover)

summaryTab = summary(InteractSigBPover)
summaryTab[summaryTab$P <= 0.05, ]

#get("GO:0006550", GOTERM)

unionRes_set1 <- summaryTab # union_GO_subset1.txt
write.table(unionRes_set1, file = "union_GO_subset1.txt", sep = '\t')

a <- unionSet_2

probeIds <- rownames(rat.rma)
GOannot = mget(probeIds, hgu133plus2GO)
notEmpty = function(x) {
	if (length(x) == 1 && is.na(x))
		FALSE
	else TRUE
}

haveGo = sapply(GOannot,notEmpty)
GoFeatures = probeIds[haveGo]

entrezIds <- mget(GoFeatures, hgu133plus2ENTREZID)
haveEntrezId <- sapply(entrezIds, function(x) !is.na(x))
entrezFeatures = GoFeatures[haveEntrezId]

require(genefilter)
EntrezSelect = names(findLargest(a, rnorm(length(a)), "hgu133plus2.db"))
EntrezUniverse = names(findLargest(entrezFeatures, rnorm(length(entrezFeatures)), "hgu133plus2.db"))


InteractBPover = new("GOHyperGParams", geneIds = EntrezSelect, universeGeneIds = EntrezUniverse, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
InteractSigBPover = hyperGTest(InteractBPover)

summaryTab = summary(InteractSigBPover)
summaryTab[summaryTab$P <= 0.05, ]

#get("GO:0006550", GOTERM)

unionRes_set2 <- summaryTab # union_GO_subset2.txt
write.table(unionRes_set2, file = "union_GO_subset2.txt", sep = '\t')


#After heatmap with interSet

interSet_1 <- c("210271_at", "203576_at", "1562828_at", "212288_at", "200845_s_at", "201416_at", "212713_at", "1555392_at", "206456_at", "205667_at", "203607_at", "200807_s_at", "204546_at", "206186_at", "202260_s_at")

interSet_2 <- c("208692_at", "204105_s_at", "202719_s_at", "200715_x_at", "206017_at", "1559149_at", "1561341_at", "205637_s_at", "212774_at", "205493_s_at", "201209_at", "200903_s_at", "205378_s_at", "204457_s_at")

a <- interSet_1

probeIds <- rownames(rat.rma)
GOannot = mget(probeIds, hgu133plus2GO)
notEmpty = function(x) {
	if (length(x) == 1 && is.na(x))
		FALSE
	else TRUE
}

haveGo = sapply(GOannot,notEmpty)
GoFeatures = probeIds[haveGo]

entrezIds <- mget(GoFeatures, hgu133plus2ENTREZID)
haveEntrezId <- sapply(entrezIds, function(x) !is.na(x))
entrezFeatures = GoFeatures[haveEntrezId]

require(genefilter)
EntrezSelect = names(findLargest(a, rnorm(length(a)), "hgu133plus2.db"))
EntrezUniverse = names(findLargest(entrezFeatures, rnorm(length(entrezFeatures)), "hgu133plus2.db"))


InteractBPover = new("GOHyperGParams", geneIds = EntrezSelect, universeGeneIds = EntrezUniverse, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
InteractSigBPover = hyperGTest(InteractBPover)

summaryTab = summary(InteractSigBPover)
summaryTab[summaryTab$P <= 0.05, ]

#get("GO:0006550", GOTERM)

interRes_set1 <- summaryTab # union_GO_subset1.txt
write.table(interRes_set1, file = "inter_GO_subset1.txt", sep = '\t')

a <- interSet_2

probeIds <- rownames(rat.rma)
GOannot = mget(probeIds, hgu133plus2GO)
notEmpty = function(x) {
	if (length(x) == 1 && is.na(x))
		FALSE
	else TRUE
}

haveGo = sapply(GOannot,notEmpty)
GoFeatures = probeIds[haveGo]

entrezIds <- mget(GoFeatures, hgu133plus2ENTREZID)
haveEntrezId <- sapply(entrezIds, function(x) !is.na(x))
entrezFeatures = GoFeatures[haveEntrezId]

require(genefilter)
EntrezSelect = names(findLargest(a, rnorm(length(a)), "hgu133plus2.db"))
EntrezUniverse = names(findLargest(entrezFeatures, rnorm(length(entrezFeatures)), "hgu133plus2.db"))


InteractBPover = new("GOHyperGParams", geneIds = EntrezSelect, universeGeneIds = EntrezUniverse, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
InteractSigBPover = hyperGTest(InteractBPover)

summaryTab = summary(InteractSigBPover)
summaryTab[summaryTab$P <= 0.05, ]

#get("GO:0006550", GOTERM)

interRes_set2 <- summaryTab # union_GO_subset2.txt
write.table(interRes_set2, file = "inter_GO_subset2.txt", sep = '\t')
```



# Annotation (my initial attempt)
```{r step 10 Gene Ontology, warning=FALSE, message=FALSE}
best2 <- read.csv("bestGenes.csv")

map <- getAnnMap("ENTREZID", "hgu133plus2", load=TRUE, type=c("env", "db"))

ll <- getEG(as.character(best2[, 1]), "hgu133plus2.db")
sym <- getSYMBOL(as.character(best2[, 1]), "hgu133plus2.db")

tab <- data.frame(sym)

# tab2 <- data.frame(sym, signif(tab[, -1], 3))

tab <- data.frame(rownames(tab), tab)
colnames(tab)[1] <- c("Probe ID") 

tab <- tab[!is.na(tab)]


```


## **Session Info** ##
```{r, warning=FALSE, message=FALSE, error=FALSE}
sessionInfo()

gc()
```
