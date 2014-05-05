---
title: "Project_Hannier&Yong_Github_042514"
author: "Hannier Pulido"
date: "April 25, 2014"
output: pdf_document
---

# **Objective**
The aim of this project is to identify a gene signature for the four different types of medulloblastoma described by (Robinson et al., 2012). The analysis reveals a list of genes whose expression is correlated with the medulloblastoma subgroups

# **Introduction**
In order to classify the samples in the four different types of medulloblastoma, we used the expression profiling by array available at the Gene Expression Omnibus website with the accession number GSE37418. A total of 54,675 genes from 76 samples were analyzed using Bioconductor and Random Forest (Breiman, 2001) for classification. 

(Develop more introduction with the results obtained from the annotation and gene ontology… what did we discover? What is the gene signature for each of the four groups?)

```{r}
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

## **2.1 Boxplots of raw data**
```{r }
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

D <- data.frame(sampleNames(CELdata))
D$genLabel <- lapply(strsplit(as.character(D$sampleNames.CELdata.), "\\_"), "[", 1) 

par(mar = c(5.1, 4.1, 4.1, 2.1))
boxplot(CELdata, col = color, las = 3, xaxt = "n", main = "Boxplot of raw data")
axis(1, at = 1:76, labels = D$genLabel, las = 3, cex.axis = 0.6)
legend(-2, 16, horiz = TRUE, c("G3", "G4", "SHH", "SHH_OUTLIER", "U", "WNT"), 
       fill = c("green3", "cyan2", "orange", "yellow","magenta", "firebrick1") )
```

## **2.2 Density plots of raw data**
```{r }
hist(CELdata, type = "l", col = color, main = "Histogram of raw data")
legend("topright", c("G3", "G4", "SHH", "SHH_OUTLIER", "U", "WNT"), 
       fill = c("green3", "cyan2", "orange", "yellow","magenta", "firebrick1") )
```

## **2.3 RNA degradation Plot**
In this plot, we are looking if any arrays are really different from the others. 
```{r}
RNA.deg <- AffyRNAdeg(CELdata)
plotAffyRNAdeg(RNA.deg, col = color)
legend("topleft", c("G3", "G4", "SHH", "SHH_OUTLIER", "U", "WNT"), 
       fill = c("green3", "cyan2", "orange", "yellow","magenta", "firebrick1") )
```

# **Step 3. Preprocess**
```{r}
cancer.rma <- rma(CELdata)
dim(cancer.rma)

eset.mas5 <- mas5calls(CELdata)

# Creates a function to remove the absent genes
mascallsfilter <- function(cutoff = "A", number){
          function(x){
          sum((x) == (cutoff)) != number
          }
    } 

m5cfun <- mascallsfilter(number = dim(eset.mas5)[2])     
mcfun <- filterfun(m5cfun)
mcsub <- genefilter(eset.mas5, mcfun)
sum(mcsub) # 47330
 
cancer.sub1 <- cancer.rma[mcsub,]
dim(cancer.sub1)        # 47330    76
cancer.sub <- get.array.subset.affybatch(cancer.sub1, "subgroup", c("G3", "G4", "SHH", "WNT"))
dim(cancer.sub)

C <- data.frame(sampleNames(cancer.sub))
C$genLabel <- lapply(strsplit(as.character(C$sampleNames.cancer.sub.), "\\_"), "[", 1) 

sampleNames(cancer.sub) <- as.character(C$genLabel)
```

# **Step 4. Preliminary**
```{r}
pdata <- pData(cancer.sub)
expr <- data.frame(exprs(cancer.sub)) # Expression data
label <- as.vector(pData(cancer.sub)$subgroup) # character vector for the subgroups
cluster_2 <- kmeans(t(expr), iter.max = 1000, 2)$cluster
res_2 <- cbind(pdata, cluster_2)
table(res_2$subgroup, res_2$cluster)
cluster_4 <- kmeans(t(expr), iter.max = 1000, 4)$cluster

res_4 <- cbind(pdata, cluster_4)
table(res_4$subgroup, res_4$cluster)

fisher.test(table(res_2$subgroup, res_2$cluster_2), workspace=2e9)
fisher.test(table(res_4$subgroup, res_4$cluster_4), workspace=2e9)

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

dist_res <- dist(t(expr))
mds<-cmdscale(dist_res) # Principal coordinates analysis for the distance matrix
x<-mds[,1]
y<-mds[,2]

plot(x, y, type = 'n', main = "MDS plot according to a k-means result")
text(x, y, label, col = color, cex = 1, font = 2)
```

## **5.1 Boxplots of filtered data**
```{r}
boxplot(exprs(cancer.sub), col = color, las = 3, xaxt = "n", main = "Boxplot of filtered data")
axis(1, at = 1:73, labels = C$genLabel, las = 3, cex.axis = 0.6)
legend(20, 15, horiz = TRUE, bg = "white", c("G3", "G4", "SHH", "WNT"), 
       fill = c("green3", "cyan2", "orange", "firebrick1") )
```

## **5.2 Density plots of filtered data**
```{r}
plot(density(expr[, 1]), col = color[[1]], ylim = c(0, 0.25), xlab = "", 
     main = "Histogram of filtered data")
for(i in 2:ncol(expr)) lines(density(expr[, i]), col = color[i])
legend("topright", c("G3", "G4", "SHH", "WNT"), 
       fill = c("green3", "cyan2", "orange", "firebrick1") )
```

## **5.3 MVA plots for filtered data**
```{r}
which(pdata$subgroup == "G3") # Returns the number of columns with samples G3
which(pdata$subgroup == "G4") # Returns the number of columns with samples G4
which(pdata$subgroup == "SHH") # Returns the number of columns with samples SHH
which(pdata$subgroup == "WNT") # Returns the number of columns with samples WNT

mva.pairs(expr[, c(9, 23, 46, 61)], col = color, cex = 1, 
          main = "MVA plot for 4 arrays of G3") 

mva.pairs(exprs(cancer.sub)[, c(1, 33, 27, 56)], col = color, cex = 1, 
          main = "MVA plot for 4 arrays of G4") 

mva.pairs(exprs(cancer.sub)[, c(5, 31, 42, 71)], col = color, cex = 1, 
          main = "MVA plot for 4 arrays of SHH") 

mva.pairs(exprs(cancer.sub)[, c(3, 47, 62, 64)], col = color, cex = 1,
          main = "MVA plot for 4 arrays of WNT") 
```

## **6.1 Hierarchical clustering**
```{r }
A <- t(expr)
dim(A)

label <- as.vector(pData(cancer.sub)$subgroup)

my.dist <- function(x) dist(x, method = "euclidean")

my.hclust <- function(d) hclust(d, method = "ward")

d <- as.dist(dist(A, method = "euclidean"))
fit.h <- hclust(d, method = "ward")
plotColoredClusters(fit.h, labs = label, cols = color, cex = 0.8, 
                    xlab = "", sub = "", main = "Hierarchical clustering")
```

## **6.2 PCA**
```{r}
pca <- prcomp(A, scale = T)
summary(pca)
scatterplot3d(pca$x[, 1:3], pch = 20, cex.symbols = 2, color = color, 
              main = "PCA for the filtered genes")
legend(5.5, 3, c("G3", "G4", "SHH", "WNT"), 
       fill = c("green3", "cyan2", "orange", "firebrick1") )
```



# **Step 7. Random Forest**
```{r}
subgroup <- factor(pdata$subgroup)
rf <- randomForest(t(exprs(cancer.sub)), factor(subgroup), mtry = 217, ntree = 500,
                   importance = TRUE, proximity = TRUE) 

colRF <- c("black", "green3", "cyan2", "orange", "firebrick1")
plot(rf, col = colRF, main = "random Forest error rate")
legend("topright", colnames(rf$err.rate), col = colRF, fill = colRF)

rf.vs1 <- varSelRF(t(exprs(cancer.sub)), factor(subgroup))

# VarSelRF
Best.vars <- data.frame(rf.vs1$selected.vars)
colnames(Best.vars) <- "Best"
head(Best.vars)
dim(Best.vars) # 58 genes are selected by VarSelRF


# Boruta 
Bor <- Boruta(t(expr), factor(subgroup), pValue = 0.05)

# Extract names of the selected attributes
Best.Bor <- data.frame(getSelectedAttributes(Bor, withTentative = FALSE))
colnames(Best.Bor) <- "Best"
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
```{r}
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
```{r}
bestGenes <- read.csv("bestGenes.csv")
B <- exprs(cancer.sub)
dim(B)

keep <- subset(B, rownames(B) %in% bestGenes[,1])
dim(keep)

bt <- t(keep)

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
```{r}
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

# Construction of a gene list that we want to enrich into GO
EntrezSelect = names(findLargest(a, rnorm(length(a)), "hgu133plus2.db"))

# Construction of the entire gene list in the microarray platform
EntrezUniverse = names(findLargest(entrezFeatures, rnorm(length(entrezFeatures)), "hgu133plus2.db"))


InteractBPover = new("GOHyperGParams", geneIds = EntrezSelect, universeGeneIds = EntrezUniverse, annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
InteractSigBPover = hyperGTest(InteractBPover)

summaryTab = summary(InteractSigBPover)
summaryTab[summaryTab$P <= 0.05, ]

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

unionRes <- summaryTab # union_GO.txt
write.table(unionRes, file = "union_GO.txt", sep = '\t')

#After heatmap with unionSet

# Subset 1 of union gene set
unionSet_1 <- c("206456_at", "218736_s_at", "220111_s_at", "226796_at", "241505_at", "218479_s_at", "224539_s_at", "242527_at", "215047_at", "215531_s_at", "241881_at", "222725_s_at", "232069_at", "235044_at", "228665_at", "212713_at", "243428_at", "203576_at", "205667_at", "243435_at", "215028_at", "204457_s_at", "232286_at", "241411_at", "201209_at", "200903_s_at", "221923_s_at", "230765_at", "1559149_at", "235004_at", "206186_at", "242959_at", "227441_s_at", "227290_at")

# Subset 1 of union gene set
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

unionRes_set2 <- summaryTab # union_GO_subset2.txt
write.table(unionRes_set2, file = "union_GO_subset2.txt", sep = '\t')

#After heatmap with interSet

# Subset 1 of intersection gene set
interSet_1 <- c("210271_at", "203576_at", "1562828_at", "212288_at", "200845_s_at", "201416_at", "212713_at", "1555392_at", "206456_at", "205667_at", "203607_at", "200807_s_at", "204546_at", "206186_at", "202260_s_at")

# Subset 2 of intersection gene set
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

interRes_set2 <- summaryTab # union_GO_subset2.txt
write.table(interRes_set2, file = "inter_GO_subset2.txt", sep = '\t')
```


## **Session Info** ##
```{r}
sessionInfo()

gc()
```
