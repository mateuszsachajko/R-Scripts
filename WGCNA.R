library(WGCNA)
library(DESeq2)
library(plyr)
library(gplots)

options(stringsAsFactors = T)
#### DATA PREP #######

datExpr = read.csv("/media/mateusz/daneMateusz/miRNASwinie/fasta/sscMir.csv")
  # head(datExpr)
  # datExpr[1,] #first row
# datExpr[,1] #first col
datTraits = read.csv("/media/mateusz/daneMateusz/miRNASwinie/phenotype/phenotypeDesc.csv")
rownames(datTraits)=datTraits[,1]
datTraits=datTraits[,c(-1,-2,-3,-4)]

datExpr=datExpr[,c(-2,-3,-4,-17:-28)]
# dim(datExpr)

#getting rid off mirs with 0 reads

keep = rowSums(datExpr[,-1])>0

datExpr = datExpr[keep,]
datExpr
table(duplicated(datExpr[,1]))


#merging data of expressed miRNA from two diferent precursors
datExpr=ddply(datExpr,"X.miRNA",numcolwise(sum))
rownames(datExpr)=datExpr[,1]
datExpr=datExpr[,-1]
datExpr[c(7,8),]

names(datExpr) = c("1R1_PxD_ctrl", "2R1_PxD_treat","3R1_P_ctrl","4R1_PxD_treat","1R2_PxD_ctrl","3R2_PxD_ctrl","30R2_PxD_ctrl",
                   "21R2_PxD_treat","32R3_P_ctrl","4R3_PxD_treat","20R3_P_treat","9R3_P_treat")
datExpr
#filtering data with low counts
keep = rowSums(datExpr>0)>=6
datExpr = datExpr[keep,]
datExpr=t(datExpr)
datExpr1 = datExpr + 1
# datExpr1 = datExpr1[,-95]

datExpr = datExpr1
# rm(datExpr1)

#normalizations

### should be varianceStabilizingTransformation(as.matrix(datExpr), blind = T, fitType = "local")
### but WGCNA still applies this type of fit automatically for this data (warning message), blind = T is set as default.
datExpr=varianceStabilizingTransformation(as.matrix(datExpr), blind = T, fitType =  "local")
# dE = log2(datExpr1)
# log2(2+1)
# #good samples genes
datExpr0=datExpr
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
  {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
# sample Tree


datTraits = cbind(GROUP = c(1,2,3,2,1,1,1,2,3,2,4,4),BREED = c(1,1,2,1,1,1,2,1,2,1,2,2),
                  DIET = c(1,2,1,2,1,1,1,2,1,2,2,2), datTraits)

# datTraits = datTraits[,-1]

sampleTree = hclust(dist(datExpr))
traitColors = numbers2colors(datTraits, signed = FALSE);
png("plots/SampleDendrotraitHM.png", width = 1920, height = 1920)
plotDendroAndColors(sampleTree, traitColors,groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",cex.rowText = 2,cex.main = 4,
                    cex.colorLabels = 2,cex.dendroLabels = 2, cex.lab = 2, cex.axis = 2, mar = c(5,12,4,2))
dev.off()
###### WGCNA  ######
getwd()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 3,blockSize = 12000)
# Plot the results:
#sizeGrWindow(9, 5)
png("plots/sftThreshold.png",width = 1920,height = 1920)
par(mfrow = c(1,2), mar =c(6,4,4,2),cex=4);
cex1 = 1;
# Scale-free topology fit index as a function of the soft-thresholding power

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),ylim=c(-0.2,1));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#options(stringsAsFactors = FALSE)
dev.off()
net = blockwiseModules(datExpr, power = 18, networkType = "signed" , corType = "bicor",
                       TOMType = "signed", labelUnlabeled = TRUE,
                       numericLabels = FALSE, pamRespectsDendro = TRUE,
                       loadTOM = F, saveTOMs = FALSE, 
                       checkMissingData = T, verbose = 3, robustY=FALSE, mergeCutHeight = 0.15,minModuleSize = 3)

#load(file = "net.RData")
list = ls()
save(list, file = "net.Rdata")

dynamicColors = labels2colors(net$colors)
#dynamicColors = dynamicColors[net$blockGenes[[(1)]]]
#print("gene tree")
geneTree = net$dendrograms[[1]]


# 
# pdf("gene_clusteringTOM.pdf",height=8,width=12) 
# sizeGrWindow(12,9)
# plot(geneTree, xlab="", sub="", 
#      main = "Gene clustering on TOM-based dissimilarity",
#      labels = FALSE, hang = 0.04);
# dev.off()
# Plot the dendrogram and the module colors underneath

png("plots/gene_clust.png",width=1920,height=1080) 
par(cex=5);
# sizeGrWindow(12,9)
plotDendroAndColors(geneTree, dynamicColors[net$blockGenes[[(1)]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, cex.main = 4,cex.lab = 2, cex.axis = 2, cex.sub = 2)
dev.off()

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")

nGenes= ncol(datExpr)
nSamples = nrow(datExpr)
# ??orderMEs
# MEs1 = orderMEs(MEs)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1), ")", sep = "")
png("plots/labeledHM.png",width=1920,height=1080) 
# sizeGrWindow(12,12)
par(cex = 1, mar = c(5, 4, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(datTraits),yLabels = names(MEs),
               colorLabels = FALSE,colors = blueWhiteRed(50),
               textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 1,zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

######## TUTAJ SKONCZONE, wyszlo spoko w miare!!!!


####### Heatmap ######

table(net$colors)

### colOrder wziac z dynamicColors
colOrder = sort(net$color)
colOrder = sort(dynamicColors)

ordergenes = order(dynamicColors)
hmDatExpr=t(datExpr)
hmDatExpr = data.frame(hmDatExpr)
scaleData = scale(log(hmDatExpr))
names(hmDatExpr) = rownames(datExpr)

png("plots/hmAllGenes.png", width = 1920, height = 1080)
heatmap.2(as.matrix(hmDatExpr), labRow = NA, margins = c(22,10),
          cexCol = 3, cexRow = 0.5, scale ="row", RowSideColors = colOrder,
          col=greenred(100), trace = "none", density.info = "none", dendrogram='row')
dev.off()

#### all sig modules
#### black, blue, brown, green, magenta, pink, purple, red, turquoise, yellow
#### all modules are significant

#### gr sl lop

keep = dynamicColors=='black'|dynamicColors=='blue'|dynamicColors=='brown'|dynamicColors=='magenta'|
dynamicColors=='purple'|dynamicColors=='turquoise'|dynamicColors=='yellow'
keep
hmDatExpr1 = hmDatExpr[keep,]
png("plots/hmGrSlLop.png", width = 1920, height = 1080)
heatmap.2(as.matrix(hmDatExpr1), labRow = NA, margins = c(22,5),
          cexCol = 3, cexRow = 0.5, scale ="row", RowSideColors = dynamicColors[keep],
          col=greenred(100), trace = "none", density.info = "none", dendrogram='row')
dev.off()

#### PE24

keep = dynamicColors=='pink'|dynamicColors=='blue'
keep
hmDatExpr1 = hmDatExpr[keep,]
png("plots/hmPE24.png", width = 1920, height = 1080)
heatmap.2(as.matrix(hmDatExpr1), labRow = NA, margins = c(22,5),
          cexCol = 3, cexRow = 0.5, scale ="row", RowSideColors = dynamicColors[keep],
          col=greenred(100), trace = "none", density.info = "none", dendrogram='row')
dev.off()

#### a*

keep = dynamicColors=='black'|dynamicColors=='blue'|dynamicColors=='brown'|dynamicColors=='magenta'|
dynamicColors=='purple'|dynamicColors=='turquoise'|dynamicColors=='yellow'|dynamicColors=='red'|dynamicColors=='green'
hmDatExpr1 = hmDatExpr[keep,]
png("plots/hmAasterix.png", width = 1920, height = 1080)
heatmap.2(as.matrix(hmDatExpr1), labRow = NA, margins = c(22,5),
          cexCol = 3, cexRow = 0.5, scale ="row", RowSideColors = dynamicColors[keep],
          col=greenred(100), trace = "none", density.info = "none", dendrogram='row')
dev.off()


#### popiol

keep = dynamicColors=='blue'|dynamicColors=='magenta'|
  dynamicColors=='purple'|dynamicColors=='turquoise'|dynamicColors=='yellow'|dynamicColors=='red'|
  dynamicColors=='green'
hmDatExpr1 = hmDatExpr[keep,]
png("plots/popiol.png", width = 1920, height = 1080)
heatmap.2(as.matrix(hmDatExpr1), labRow = NA, margins = c(22,5),
          cexCol = 3, cexRow = 0.5, scale ="row", RowSideColors = dynamicColors[keep],
          col=greenred(100), trace = "none", density.info = "none", dendrogram='row')
dev.off()

#### C20.5n3
  
keep = dynamicColors=='red'|dynamicColors=='green'
keep
hmDatExpr1 = hmDatExpr[keep,]
png("plots/hmC20_5n3.png", width = 1920, height = 1080)
heatmap.2(as.matrix(hmDatExpr1), labRow = NA, margins = c(22,5),
          cexCol = 3, cexRow = 0.5, scale ="row", RowSideColors = dynamicColors[keep],
          col=greenred(100), trace = "none", density.info = "none", dendrogram='row')
dev.off()
