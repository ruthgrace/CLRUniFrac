#!/usr/bin/env Rscript

# consider testing with different rarification methods later ?

options(error=recover)

library(ape)
library(phangorn)
library(vegan)

originalPar <- par()


source("../metrics.r")

replicates <- 5
extraReplicates <- 10

darkorchid <- col2rgb("darkorchid4")
transparentdarkorchid <- rgb(darkorchid[1]/255,darkorchid[2]/255,darkorchid[3]/255,0.3)

aquamarine <- col2rgb("aquamarine4")
transparentaquamarine <- rgb(aquamarine[1]/255,aquamarine[2]/255,aquamarine[3]/255,0.3)

pcoaLabels <- c(rep("pcoa1",3),rep("pcoa12",3),rep("pcoa123",3))
unifracLabels <- rep(c("uwUnifrac","wUnifrac","iUnifrac"),3)
plotDataColNames <- paste(pcoaLabels, unifracLabels, sep = ".")

pcoaText <- c(rep("PCoA first component",3),rep("PCoA first 2 components",3),rep("PCoA first 3 components",3))
unifracText <- rep(c("Unweighted UniFrac","Weighted UniFrac","Information UniFrac"),3)
plotDataTextLabels <- paste(unifracText,pcoaText, sep = "\n")

screeText <- paste("Axis",rep(c(1:5),3))
screeUnifracTest <- rep(c("Unweighted UniFrac","Weighted UniFrac","Information UniFrac"),5)
screeLabels <- paste(screeUnifracTest,screeText)



rootTree <- function(tree) {
	if (!is.rooted(tree)) {
		tree <- midpoint(tree)
	}
	return(tree)
}

removeTreeTipLabelSingleQuotes <- function(tree) {
	tree$tip.label <- gsub("'","",tree$tip.label)
	return(tree)
}

runReplicate <- function(otu,groups,tree,nSamples) {
	
	
	#sample 50 samples from condition 1
	group1.indices <- which(groups==levels(groups)[1])
	group1.rand <- otu[as.integer(sample(group1.indices,nSamples,replace=FALSE)),]
	
	#sample 50 samples from condition 2
	group2.indices <- which(groups==levels(groups)[2])
	group2.rand <- otu[as.integer(sample(group2.indices,nSamples,replace=FALSE)),]
	
	#concatenate
	data <- rbind(group1.rand,group2.rand)

	#make groups
	newGroups <- c(rep(levels(groups)[1],nSamples),rep(levels(groups)[2],nSamples))
	return(getAllPcoaMetrics(data,newGroups,tree))
}

runMixedReplicate <- function(otu1,otu2,groups1,groups2,tree,nSamples) {	
	#sample 50 samples from condition 1,group1
	group1.indices <- which(groups1==levels(groups1)[1])
	group1.rand <- otu1[as.integer(sample(group1.indices,nSamples,replace=FALSE)),]
	
	#sample 50 samples from condition 2,group2
	group2.indices <- which(groups2==levels(groups2)[2])
	group2.rand <- otu2[as.integer(sample(group2.indices,nSamples,replace=FALSE)),]
	
	#concatenate
	data <- rbind(group1.rand,group2.rand)

	#make groups
	newGroups <- c(rep(levels(groups1)[1],nSamples),rep(levels(groups2)[2],nSamples))
	return(getAllPcoaMetrics(data,newGroups,tree))
}

addDisimilarOTUs <- function(otu1,otu2,tree) {
	# remove OTUs not in tree
	otu1.tree.otu <- otu1[,which(colnames(otu1) %in% tree$tip.label)]
	otu2.tree.otu <- otu2[,which(colnames(otu2) %in% tree$tip.label)]

	# find indices of otus common and unique to both lists
	otu1.common.otus <- which(colnames(otu1.tree.otu) %in% colnames(otu2.tree.otu))
	otu2.common.otus <- which(colnames(otu2.tree.otu) %in% colnames(otu1.tree.otu))
	otu1.unique.otus <- which(!(colnames(otu1.tree.otu) %in% colnames(otu1.tree.otu)[otu1.common.otus]))
	otu2.unique.otus <- which(!(colnames(otu2.tree.otu) %in% colnames(otu2.tree.otu)[otu2.common.otus]))

	# construct blank otus from otu2 to add to otu1
	otu2.unique.blanks <- otu2.tree.otu[,otu2.unique.otus]
	otu2.unique.blanks[,] <- 0
	if (nrow(otu2)>nrow(otu1)) {
		otu2.unique.blanks <- otu2.unique.blanks[1:nrow(otu1),]
	}
	else {
		otu2.unique.blanks <- otu2.unique.blanks[c(c(1:nrow(otu2)),rep(1,(nrow(otu1)-nrow(otu2)))),]
	}
	newotu1 <- data.frame(otu1.tree.otu,otu2.unique.blanks)

	# construct blank otus from otu1 to add to otu2
	otu1.unique.blanks <- otu1.tree.otu[,otu1.unique.otus]
	otu1.unique.blanks[,] <- 0
	if (nrow(otu1)>nrow(otu2)) {
		otu1.unique.blanks <- otu1.unique.blanks[1:nrow(otu2),]
	}
	else {
		otu1.unique.blanks <- otu1.unique.blanks[c(c(1:nrow(otu1)),rep(1,(nrow(otu2)-nrow(otu1)))),]
	}
	newotu2 <- data.frame(otu2.tree.otu,otu1.unique.blanks[1:nrow(otu2),])

	#return both new otu tables
	returnList <- list()
	returnList[[1]] <- newotu1
	returnList[[2]] <- newotu2
	return(returnList)
}

compare <- function(replicateMethod, otuList, groupList, comparisonList, treeList, comparisonTitleList, fileName, nSamples) {
	pdf(fileName)
	print(paste("length of comparisonList",length(comparisonList)))
	for (i in 1:length(comparisonList)) {
		print(paste("comparison list",i))
		print(comparisonList[[i]])
		compare <- comparisonList[[i]]
		index1 <- compare[1]
		otu1 <- otuList[[index1]]
		group1 <- groupList[[index1]]
		plotTitle <- comparisonTitleList[i]
		tree <- treeList[[i]]
		if (identical(replicateMethod,runReplicate)) {
			getReplicate(replicateMethod,otu1,group1,tree,plotTitle,nSamples)
		}
		else {
			index2 <- compare[2]
			otu2 <- otuList[[index2]]
			group2 <- groupList[[index2]]
			getReplicate(replicateMethod,otu1,group1,tree,plotTitle,nSamples,group2=group2,otu2=otu2)
		}
	}

	dev.off()

}

getReplicate <- function(replicateMethod,otu1,group1,tree,plotTitle,nSamples,group2=NULL,otu2=NULL) {
	#do comparisons
	reps <- list()
	#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
	plot.data <- data.frame(matrix(nrow=5,ncol=9))
	colnames(plot.data) <- plotDataTextLabels

	dist <- data.frame(matrix(nrow=5,ncol=9))
	colnames(dist) <- plotDataTextLabels
	sd <- data.frame(matrix(nrow=5,ncol=9))
	colnames(sd) <- plotDataTextLabels

	scree <- data.frame(matrix(nrow=5,ncol=15))
	colnames(scree) <- screeLabels


	for (i in 1:replicates) {
		if (identical(replicateMethod,runReplicate)) {
			reps[[i]] <- runReplicate(otu1,group1,tree,nSamples)
		}
		else {
			reps[[i]] <- runMixedReplicate(otu1,otu2,group1,group2,tree,nSamples)
		}
		plot.data[i,] <- unlist(data.frame(t(reps[[i]]$effect)))
		dist[i,] <- c(reps[[i]]$meanDist.SD[[1]]$meanDist,reps[[i]]$meanDist.SD[[2]]$meanDist,reps[[i]]$meanDist.SD[[3]]$meanDist)
		sd[i,] <- c(reps[[i]]$meanDist.SD[[1]]$error,reps[[i]]$meanDist.SD[[2]]$error,reps[[i]]$meanDist.SD[[3]]$error)
		scree[i,] <- c(reps[[i]]$screeData$uwUnifrac[1:5],reps[[i]]$screeData$wUnifrac[1:5],reps[[i]]$screeData$eUnifrac[1:5])
	}

	#make plots
	originalPar <- par()

	#effect size plot
	par(mar=c(13, 4, 4, 2) + 0.1)
	par(cex.lab=1.3)
	par(cex.main=1.5)
	stripchart(plot.data,vertical=TRUE,main="Sparsity filter at 0.1%",group.names=colnames(plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
	par(originalPar)

	#mean separation with standard deviation plot (averaged over 5 replicates .....)
	par(mar=c(13, 6, 4, 2) + 0.1)

	meanDist <- unlist(lapply(dist,mean))
	meanSd <- unlist(lapply(sd,mean))
	myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sparsity filter at 0.1%")
	segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
	segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
	segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

	# scree plots
	par(originalPar)
	par(mar=c(12, 6, 4, 2) + 0.1)
	screePlotData <- apply(scree,2,mean)
	screePlotError <- apply(scree,2,sd)
	myBarPlot <- barplot(screePlotData,col=transparentdarkorchid,las=2,ylim=c(0,1),ylab="Variation explained by each axis of PCoA",main="Sparsity filter at 0.1%")
	segments(myBarPlot, screePlotData - screePlotError, myBarPlot, screePlotData + screePlotError, lwd=2)
	segments(myBarPlot - 0.1, screePlotData - screePlotError, myBarPlot + 0.1, screePlotData - screePlotError, lwd=2)
	segments(myBarPlot - 0.1, screePlotData + screePlotError, myBarPlot + 0.1, screePlotData + screePlotError, lwd=2)

	par(originalPar)
	# pcoa plots, only plot first replicate in each comparison
	palette(c(transparentdarkorchid,transparentaquamarine,"blue","black"))
	pcoaGroups <- as.factor(c(rep(1,nSamples),rep(2,nSamples)))
	plot(reps[[1]]$pcoa$uwUnifrac$vectors[,1],reps[[1]]$pcoa$uwUnifrac$vectors[,2], type="p",col=pcoaGroups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(reps[[1]]$screeData$uwUnifrac[1],digits=3),"variance explained"),ylab=paste("Second Component", round(reps[[1]]$screeData$uwUnifrac[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(reps[[1]]$pcoa$wUnifrac$vectors[,1],reps[[1]]$pcoa$uwUnifrac$vectors[,2], type="p",col=pcoaGroups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(reps[[1]]$screeData$wUnifrac[1],digits=3),"variance explained"),ylab=paste("Second Component", round(reps[[1]]$screeData$wUnifrac[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(reps[[1]]$pcoa$eUnifrac$vectors[,1],reps[[1]]$pcoa$uwUnifrac$vectors[,2], type="p",col=pcoaGroups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(reps[[1]]$screeData$eUnifrac[1],digits=3),"variance explained"),ylab=paste("Second Component", round(reps[[1]]$screeData$eUnifrac[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)


	par(originalPar)
	# return replicate
	return(reps)
}

#all CLR DIRICHLET commented out while the method is being fixed.
#attempt at sparsity filter instead -- 30 counts minimum is what is stable


low.data <- read.table("low_sequencing_depth_hmp_data.txt",sep="\t",header=TRUE,row.names=1)
med.data <- read.table("med_sequencing_depth_hmp_data.txt",sep="\t",header=TRUE,row.names=1)
high.data <- read.table("high_sequencing_depth_hmp_data.txt",sep="\t",header=TRUE,row.names=1)

#get metadata (stool vs saliva)
low.groups <- as.factor(gsub("_.*", "", colnames(low.data)))
med.groups <- as.factor(gsub("_.*", "", colnames(med.data)))
high.groups <- as.factor(gsub("_.*", "", colnames(high.data)))

low.tree <- read.tree("./low_sequencing_depth_subtree.tre")
med.tree <- read.tree("./med_sequencing_depth_subtree.tre")
high.tree <- read.tree("./high_sequencing_depth_subtree.tre")

#get rid of extra quotes on OTU labels
low.tree <- removeTreeTipLabelSingleQuotes(low.tree)
med.tree <- removeTreeTipLabelSingleQuotes(med.tree)
high.tree <- removeTreeTipLabelSingleQuotes(high.tree)

#root tree by midpoint if not rooted
low.tree <- rootTree(low.tree)
med.tree <- rootTree(med.tree)
high.tree <- rootTree(high.tree)

#source("../../CLRUniFrac.R")
source("../../GUniFrac.R")
source("../../EntropyUniFrac.R")
#source("../../CLRDirichletUniFrac.R")

#format otu table for input into unifrac methods (rownames are samples, colnames are OTUs)
low.data.t <- t(low.data)
med.data.t <- t(med.data)
high.data.t <- t(high.data)

#get rid of any OTUs that aren't in the tree (only a couple reads discarded in total)
low.otu.unordered <- low.data.t[,which(colnames(low.data.t) %in% low.tree$tip.label)]
med.otu.unordered <- med.data.t[,which(colnames(med.data.t) %in% med.tree$tip.label)]
high.otu.unordered <- high.data.t[,which(colnames(high.data.t) %in% high.tree$tip.label)]

#order otus by abundance (least to most)
taxaOrder <- rev(order(apply(low.otu.unordered,2,sum)))
low.otu <- low.otu.unordered[,taxaOrder]
taxaOrder <- rev(order(apply(med.otu.unordered,2,sum)))
med.otu <- med.otu.unordered[,taxaOrder]
taxaOrder <- rev(order(apply(high.otu.unordered,2,sum)))
high.otu <- high.otu.unordered[,taxaOrder]



selfComparisonList3 <- list()
selfComparisonList3[[1]] <- c(1)
selfComparisonList3[[2]] <- c(2)
selfComparisonList3[[3]] <- c(3)

mixedComparisonList3 <- list()
mixedComparisonList3[[1]] <- c(1,2)
mixedComparisonList3[[2]] <- c(1,3)
mixedComparisonList3[[3]] <- c(2,3)

selfComparisonList2 <- list()
selfComparisonList2[[1]] <- c(1)
selfComparisonList2[[2]] <- c(2)

mixedComparisonList2 <- list()
mixedComparisonList2[[1]] <- c(1,2)
mixedComparisonList2[[2]] <- c(2,1)

# sparsityLabels <- c(rep("sparsity.001",3),rep("sparsity.0001",3),rep("sparsity.00001",3))
# plotSparsityDataColNames <- paste(sparsityLabels, unifracLabels, sep = ".")

# diversityLabels <- c(rep("low.diversity",3),rep("high.diversity",3))
# shortUnifracLabels <- rep(c("uwUnifrac","wUnifrac","iUnifrac"),2)
# plotDiversityDataColNames <- paste(diversityLabels, shortUnifracLabels, sep = ".")


# SEQUENCING DEPTH TEST

depth <- list()

depth$low <- low.otu
depth$med <- med.otu
depth$high <- high.otu

depth.groups <- list()
depth.groups$low <- low.groups
depth.groups$med <- med.groups
depth.groups$high <- high.groups

depth.tree <- list()
depth.tree$low <- low.tree
depth.tree$med <- med.tree
depth.tree$high <- high.tree

comparisonTitleList <- c("Sequencing depth < 3000 reads/sample","Sequencing depth 3000-6000 reads/sample","Sequencing depth > 6000 reads/sample")
fileName <- "sequencingDepthPlots.pdf"
nSamples <- 50

compare(runReplicate, depth, depth.groups, selfComparisonList3, depth.tree, comparisonTitleList, fileName, nSamples)



#low
# low.seq.depth.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# low.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.seq.depth.plot.data) <- plotDataTextLabels

# low.seq.depth.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.seq.depth.dist) <- plotDataTextLabels
# low.seq.depth.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.seq.depth.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	low.seq.depth.reps[[i]] <- runReplicate(low.otu,low.groups,low.tree,50)
# 	low.seq.depth.plot.data[i,] <- unlist(data.frame(t(low.seq.depth.reps[[i]]$effect)))
# 	low.seq.depth.dist[i,] <- c(low.seq.depth.reps[[i]]$meanDist.SD[[1]]$meanDist,low.seq.depth.reps[[i]]$meanDist.SD[[2]]$meanDist,low.seq.depth.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	low.seq.depth.sd[i,] <- c(low.seq.depth.reps[[i]]$meanDist.SD[[1]]$error,low.seq.depth.reps[[i]]$meanDist.SD[[2]]$error,low.seq.depth.reps[[i]]$meanDist.SD[[3]]$error)
# }


# #med
# med.seq.depth.reps <- list()
# med.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(med.seq.depth.plot.data) <- plotDataTextLabels

# med.seq.depth.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(med.seq.depth.dist) <- plotDataTextLabels
# med.seq.depth.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(med.seq.depth.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	med.seq.depth.reps[[i]] <- runReplicate(med.otu,med.groups,med.tree,50)
# 	med.seq.depth.plot.data[i,] <- unlist(data.frame(t(med.seq.depth.reps[[i]]$effect)))
# 	med.seq.depth.dist[i,] <- c(med.seq.depth.reps[[i]]$meanDist.SD[[1]]$meanDist,med.seq.depth.reps[[i]]$meanDist.SD[[2]]$meanDist,med.seq.depth.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	med.seq.depth.sd[i,] <- c(med.seq.depth.reps[[i]]$meanDist.SD[[1]]$error,med.seq.depth.reps[[i]]$meanDist.SD[[2]]$error,med.seq.depth.reps[[i]]$meanDist.SD[[3]]$error)
# }
# #high
# high.seq.depth.reps <- list()
# high.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(high.seq.depth.plot.data) <- plotDataTextLabels

# high.seq.depth.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(high.seq.depth.dist) <- plotDataTextLabels
# high.seq.depth.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(high.seq.depth.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	high.seq.depth.reps[[i]] <- runReplicate(high.otu,high.groups,high.tree,50)
# 	high.seq.depth.plot.data[i,] <- unlist(data.frame(t(high.seq.depth.reps[[i]]$effect)))
# 	high.seq.depth.dist[i,] <- c(high.seq.depth.reps[[i]]$meanDist.SD[[1]]$meanDist,high.seq.depth.reps[[i]]$meanDist.SD[[2]]$meanDist,high.seq.depth.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	high.seq.depth.sd[i,] <- c(high.seq.depth.reps[[i]]$meanDist.SD[[1]]$error,high.seq.depth.reps[[i]]$meanDist.SD[[2]]$error,high.seq.depth.reps[[i]]$meanDist.SD[[3]]$error)
# }

# #plot
# pdf("sequencingDepthPlots.pdf")
# par(mar=c(13, 4, 4, 2) + 0.1)
# par(cex.lab=1.3)
# par(cex.main=1.5)
# stripchart(low.seq.depth.plot.data,vertical=TRUE,main="Sequencing depth < 3000 reads/sample",group.names=colnames(low.seq.depth.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# stripchart(med.seq.depth.plot.data,vertical=TRUE,main="Sequencing depth 3000-6000 reads/sample",group.names=colnames(med.seq.depth.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# stripchart(high.seq.depth.plot.data,vertical=TRUE,main="Sequencing depth > 6000 reads/sample",group.names=colnames(high.seq.depth.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# par(originalPar)

# par(mar=c(13, 6, 4, 2) + 0.1)

# meanDist <- unlist(lapply(low.seq.depth.dist,mean))
# meanSd <- unlist(lapply(low.seq.depth.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sequencing depth < 3000 reads/sample")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# meanDist <- unlist(lapply(med.seq.depth.dist,mean))
# meanSd <- unlist(lapply(med.seq.depth.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sequencing depth 3000-6000 reads/sample")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# meanDist <- unlist(lapply(high.seq.depth.dist,mean))
# meanSd <- unlist(lapply(high.seq.depth.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sequencing depth > 6000 reads/sample")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)
# par(originalPar)


# dev.off()



# SEQUENCING DEPTH DIFFERENCE TEST

depth.diff.tree <- list()
depth.diff.tree[[1]] <- high.tree
depth.diff.tree[[2]] <- high.tree
depth.diff.tree[[3]] <- high.tree

comparisonTitleList <- c("Sequencing depth < 3000 vs. 3000-6000 reads/sample","Sequencing depth 3000-6000 vs. > 6000 reads/sample","Sequencing depth 3000-6000 vs. > 6000 reads/sample")
fileName <- "sequencingDepthDifferenceTestPlots.pdf"
nSamples <- 50

compare(runMixedReplicate, depth, depth.groups, mixedComparisonList3, depth.diff.tree, comparisonTitleList, fileName, nSamples)


# #low/med
# low.seq.depth.diff.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# low.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.seq.depth.diff.plot.data) <- plotDataTextLabels
# addedOTUs <- addDisimilarOTUs(low.otu,med.otu,high.tree)
# low.with.med.otu <- addedOTUs[[1]]
# med.with.low.otu <- addedOTUs[[2]]

# low.with.med.otu.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.with.med.otu.dist) <- plotDataTextLabels
# low.with.med.otu.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.with.med.otu.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	low.seq.depth.diff.reps[[i]] <- runMixedReplicate(low.with.med.otu,med.with.low.otu,low.groups,med.groups,high.tree,50)
# 	low.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(low.seq.depth.diff.reps[[i]]$effect)))
# 	low.seq.depth.diff.dist[i,] <- c(low.seq.depth.diff.reps[[i]]$meanDist.SD[[1]]$meanDist,low.seq.depth.diff.reps[[i]]$meanDist.SD[[2]]$meanDist,low.seq.depth.diff.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	low.seq.depth.diff.sd[i,] <- c(low.seq.depth.diff.reps[[i]]$meanDist.SD[[1]]$error,low.seq.depth.diff.reps[[i]]$meanDist.SD[[2]]$error,low.seq.depth.diff.reps[[i]]$meanDist.SD[[3]]$error)

# }
# #med/high
# med.seq.depth.diff.reps <- list()
# med.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(med.seq.depth.diff.plot.data) <- plotDataTextLabels
# addedOTUs <- addDisimilarOTUs(med.otu,high.otu,high.tree)
# med.with.high.otu <- addedOTUs[[1]]
# high.with.med.otu <- addedOTUs[[2]]

# med.with.high.otu.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(med.with.high.otu.dist) <- plotDataTextLabels
# med.with.high.otu.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(med.with.high.otu.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	med.seq.depth.diff.reps[[i]] <- runMixedReplicate(med.with.high.otu,high.with.med.otu,med.groups,high.groups,high.tree,50)
# 	med.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(med.seq.depth.diff.reps[[i]]$effect)))
# 	med.seq.depth.diff.dist[i,] <- c(med.seq.depth.diff.reps[[i]]$meanDist.SD[[1]]$meanDist,med.seq.depth.diff.reps[[i]]$meanDist.SD[[2]]$meanDist,med.seq.depth.diff.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	med.seq.depth.diff.sd[i,] <- c(med.seq.depth.diff.reps[[i]]$meanDist.SD[[1]]$error,med.seq.depth.diff.reps[[i]]$meanDist.SD[[2]]$error,med.seq.depth.diff.reps[[i]]$meanDist.SD[[3]]$error)
# }
# #low/high
# high.seq.depth.diff.reps <- list()
# high.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(high.seq.depth.diff.plot.data) <- plotDataTextLabels
# addedOTUs <- addDisimilarOTUs(low.otu,high.otu,high.tree)
# low.with.high.otu <- addedOTUs[[1]]
# high.with.low.otu <- addedOTUs[[2]]

# low.with.high.otu.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.with.high.otu.dist) <- plotDataTextLabels
# low.with.high.otu.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.with.high.otu.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	high.seq.depth.diff.reps[[i]] <- runMixedReplicate(low.with.high.otu,high.with.low.otu,low.groups,high.groups,high.tree,50)
# 	high.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(high.seq.depth.diff.reps[[i]]$effect)))
# 	high.seq.depth.diff.dist[i,] <- c(high.seq.depth.diff.reps[[i]]$meanDist.SD[[1]]$meanDist,high.seq.depth.diff.reps[[i]]$meanDist.SD[[2]]$meanDist,high.seq.depth.diff.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	high.seq.depth.diff.sd[i,] <- c(high.seq.depth.diff.reps[[i]]$meanDist.SD[[1]]$error,high.seq.depth.diff.reps[[i]]$meanDist.SD[[2]]$error,high.seq.depth.diff.reps[[i]]$meanDist.SD[[3]]$error)
# }

# pdf("sequencingDepthDifferencePlots.pdf")
# par(mar=c(13, 4, 4, 2) + 0.1)
# par(cex.lab=1.3)
# par(cex.main=1.5)
# stripchart(low.seq.depth.diff.plot.data,vertical=TRUE,main="Sequencing depth < 3000 vs. 3000-6000 reads/sample",group.names=colnames(low.seq.depth.diff.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# stripchart(med.seq.depth.diff.plot.data,vertical=TRUE,main="Sequencing depth 3000-6000 vs. > 6000 reads/sample",group.names=colnames(med.seq.depth.diff.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# stripchart(high.seq.depth.diff.plot.data,vertical=TRUE,main="Sequencing depth 3000-6000 vs. > 6000 reads/sample",group.names=colnames(high.seq.depth.diff.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# par(originalPar)

# par(mar=c(13, 6, 4, 2) + 0.1)

# meanDist <- unlist(lapply(low.seq.depth.diff.dist,mean))
# meanSd <- unlist(lapply(low.seq.depth.diff.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sequencing depth < 3000 vs. 3000-6000 reads/sample")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# meanDist <- unlist(lapply(med.seq.depth.diff.dist,mean))
# meanSd <- unlist(lapply(med.seq.depth.diff.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sequencing depth 3000-6000 vs. > 6000 reads/sample")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# meanDist <- unlist(lapply(high.seq.depth.diff.dist,mean))
# meanSd <- unlist(lapply(high.seq.depth.diff.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sequencing depth 3000-6000 vs. > 6000 reads/sample")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# par(originalPar)
# dev.off()

# # SPARSITY TEST

# #remove OTUs rarer than a thresh hold throughout all samples
# high.otu.sum <- apply(high.otu,2,sum)
# high.total.sum <- sum(high.otu)

# sparse <- list()

# sparse$otu.001 <- high.otu[,(which(high.otu.sum >= (0.001*high.total.sum)))]
# sparse$otu.0001 <- high.otu[,(which(high.otu.sum >= (0.0001*high.total.sum)))]
# sparse$otu.00001 <- high.otu[,(which(high.otu.sum >= (0.00001*high.total.sum)))]

# sparse.groups <- list()
# sparse.groups[[1]] <- sparse.groups[[2]] <- sparse.groups[[3]] <- high.groups

# sparse.tree <- list()
# sparse.tree[[1]] <- sparse.tree[[2]] <- sparse.tree[[3]] <- high.tree

# comparisonTitleList <- c("Sparsity filter at 0.1%","Sparsity filter at 0.01%","Sparsity filter at 0.001%")
# fileName <- "sparsityTestPlots.pdf"
# nSamples <- 50

# compare(runReplicate, sparse, sparse.groups, selfComparisonList3, sparse.tree, comparisonTitleList, fileName, nSamples)


# #0.1% sparsity filter
# sparse.otu.001 <- high.otu[,(which(high.otu.sum >= (0.001*high.total.sum)))]
# sparse.otu.001.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# sparse.otu.001.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.otu.001.plot.data) <- plotDataTextLabels

# sparse.otu.001.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.otu.001.dist) <- plotDataTextLabels
# sparse.otu.001.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.otu.001.sd) <- plotDataTextLabels

# sparse.otu.001.scree <- data.frame(matrix(nrow=5,ncol=15))
# colnames(sparse.otu.001.scree) <- screeLabels


# for (i in 1:replicates) {
# 	sparse.otu.001.reps[[i]] <- runReplicate(sparse.otu.001,high.groups,high.tree,50)
# 	sparse.otu.001.plot.data[i,] <- unlist(data.frame(t(sparse.otu.001.reps[[i]]$effect)))
# 	sparse.otu.001.dist[i,] <- c(sparse.otu.001.reps[[i]]$meanDist.SD[[1]]$meanDist,sparse.otu.001.reps[[i]]$meanDist.SD[[2]]$meanDist,sparse.otu.001.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	sparse.otu.001.sd[i,] <- c(sparse.otu.001.reps[[i]]$meanDist.SD[[1]]$error,sparse.otu.001.reps[[i]]$meanDist.SD[[2]]$error,sparse.otu.001.reps[[i]]$meanDist.SD[[3]]$error)
# 	sparse.otu.001.scree[i,] <- c(sparse.otu.001.reps[[i]]$screeData$uwUnifrac[1:5],sparse.otu.001.reps[[i]]$screeData$wUnifrac[1:5],sparse.otu.001.reps[[i]]$screeData$eUnifrac[1:5])
# }

# #0.01% sparsity filter
# sparse.otu.0001 <- high.otu[,(which(high.otu.sum >= (0.0001*high.total.sum)))]
# sparse.otu.0001.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# sparse.otu.0001.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.otu.0001.plot.data) <- plotDataTextLabels

# sparse.otu.0001.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.otu.0001.dist) <- plotDataTextLabels
# sparse.otu.0001.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.otu.0001.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	sparse.otu.0001.reps[[i]] <- runReplicate(sparse.otu.0001,high.groups,high.tree,50)
# 	sparse.otu.0001.plot.data[i,] <- unlist(data.frame(t(sparse.otu.0001.reps[[i]]$effect)))
# 	sparse.otu.0001.dist[i,] <- c(sparse.otu.0001.reps[[i]]$meanDist.SD[[1]]$meanDist,sparse.otu.0001.reps[[i]]$meanDist.SD[[2]]$meanDist,sparse.otu.0001.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	sparse.otu.0001.sd[i,] <- c(sparse.otu.0001.reps[[i]]$meanDist.SD[[1]]$error,sparse.otu.0001.reps[[i]]$meanDist.SD[[2]]$error,sparse.otu.0001.reps[[i]]$meanDist.SD[[3]]$error)
# }

# #0.001% sparsity filter
# sparse.otu.00001 <- high.otu[,(which(high.otu.sum >= (0.00001*high.total.sum)))]
# sparse.otu.00001.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# sparse.otu.00001.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.otu.00001.plot.data) <- plotDataTextLabels

# sparse.otu.00001.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.otu.00001.dist) <- plotDataTextLabels
# sparse.otu.00001.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.otu.00001.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	sparse.otu.00001.reps[[i]] <- runReplicate(sparse.otu.00001,high.groups,high.tree,50)
# 	sparse.otu.00001.plot.data[i,] <- unlist(data.frame(t(sparse.otu.00001.reps[[i]]$effect)))
# 	sparse.otu.00001.dist[i,] <- c(sparse.otu.00001.reps[[i]]$meanDist.SD[[1]]$meanDist,sparse.otu.00001.reps[[i]]$meanDist.SD[[2]]$meanDist,sparse.otu.00001.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	sparse.otu.00001.sd[i,] <- c(sparse.otu.00001.reps[[i]]$meanDist.SD[[1]]$error,sparse.otu.00001.reps[[i]]$meanDist.SD[[2]]$error,sparse.otu.00001.reps[[i]]$meanDist.SD[[3]]$error)
# }



# ##### NEED TO MAKE AXIS LOOK RIGHT + plot other stuff

# pdf("sparsityTestPlots.pdf")
# par(mar=c(13, 4, 4, 2) + 0.1)
# par(cex.lab=1.3)
# par(cex.main=1.5)
# stripchart(sparse.otu.001.plot.data,vertical=TRUE,main="Sparsity filter at 0.1%",group.names=colnames(sparse.otu.001.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# stripchart(sparse.otu.0001.plot.data,vertical=TRUE,main="Sparsity filter at 0.01%",group.names=colnames(sparse.otu.0001.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# stripchart(sparse.otu.00001.plot.data,vertical=TRUE,main="sparsity filter at 0.001%",group.names=colnames(sparse.otu.00001.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# par(originalPar)

# par(mar=c(13, 6, 4, 2) + 0.1)

# meanDist <- unlist(lapply(sparse.otu.001.dist,mean))
# meanSd <- unlist(lapply(sparse.otu.001.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sparsity filter at 0.1%")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# meanDist <- unlist(lapply(sparse.otu.0001.dist,mean))
# meanSd <- unlist(lapply(sparse.otu.0001.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sparsity filter at 0.01%")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# meanDist <- unlist(lapply(sparse.otu.00001.dist,mean))
# meanSd <- unlist(lapply(sparse.otu.00001.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sparsity filter at 0.001%")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# par(originalPar)

# screePlotData <- apply(sparse.otu.001.scree,2,mean)
# screePlotError <- apply(sparse.otu.001.scree,2,sd)
# myBarPlot <- barplot(screePlotData,col=transparentdarkorchid,las=2,ylim=c(0,1),ylab="Variation explained by each axis of PCoA",main="Sparsity filter at 0.1%")
# segments(myBarPlot, screePlotData - screePlotError, myBarPlot, screePlotData + screePlotError, lwd=2)
# segments(myBarPlot - 0.1, screePlotData - screePlotError, myBarPlot + 0.1, screePlotData - screePlotError, lwd=2)
# segments(myBarPlot - 0.1, screePlotData + screePlotError, myBarPlot + 0.1, screePlotData + screePlotError, lwd=2)


# dev.off()

# SPARSITY DIFFERENCE TEST

sparse.diff <- list()

sparse.diff$otu.001 <- high.otu
sparse.diff$otu.001[,(which(high.otu.sum < (0.001*high.total.sum)))] <- 0
sparse.diff$otu.0001 <- high.otu
sparse.diff$otu.0001[,(which(high.otu.sum < (0.0001*high.total.sum)))] <- 0
sparse.diff$otu.00001 <- high.otu
sparse.diff$otu.00001[,(which(high.otu.sum < (0.00001*high.total.sum)))] <- 0


sparse.diff.groups <- list()
sparse.diff.groups[[1]] <- sparse.diff.groups[[2]] <- sparse.diff.groups[[3]] <- high.groups

sparse.diff.tree <- list()
sparse.diff.tree[[1]] <- sparse.diff.tree[[2]] <- sparse.diff.tree[[3]] <- high.tree

comparisonTitleList <- c("Sparsity filter at 0.1% vs 0.01%","Sparsity filter at 0.1% vs 0.001%","Sparsity filter at 0.01% vs 0.001%")
fileName <- "sparsityDifferenceTestPlots.pdf"
nSamples <- 50

compare(runMixedReplicate, sparse.diff, sparse.diff.groups, mixedComparisonList3, sparse.diff.tree, comparisonTitleList, fileName, nSamples)



# #low/med
# sparse.diff.otu.1.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# sparse.diff.otu.1.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.diff.otu.1.plot.data) <- plotDataTextLabels

# sparse.diff.1.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.diff.1.dist) <- plotDataTextLabels
# sparse.diff.1.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.diff.1.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	sparse.diff.otu.1.reps[[i]] <- runMixedReplicate(sparse.diff.otu.001,sparse.diff.otu.0001,high.groups,high.groups,high.tree,50)
# 	sparse.diff.otu.1.plot.data[i,] <- unlist(data.frame(t(sparse.diff.otu.1.reps[[i]]$effect)))
# 	sparse.diff.1.dist[i,] <- c(sparse.diff.1.reps[[i]]$meanDist.SD[[1]]$meanDist,sparse.diff.1.reps[[i]]$meanDist.SD[[2]]$meanDist,sparse.diff.1.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	sparse.diff.1.sd[i,] <- c(sparse.diff.1.reps[[i]]$meanDist.SD[[1]]$error,sparse.diff.1.reps[[i]]$meanDist.SD[[2]]$error,sparse.diff.1.reps[[i]]$meanDist.SD[[3]]$error)

# }

# #low/high
# sparse.diff.otu.2.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# sparse.diff.otu.2.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.diff.otu.2.plot.data) <- plotDataTextLabels

# sparse.diff.2.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.diff.2.dist) <- plotDataTextLabels
# sparse.diff.2.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.diff.2.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	sparse.diff.otu.2.reps[[i]] <- runMixedReplicate(sparse.diff.otu.001,sparse.diff.otu.00001,high.groups,high.groups,high.tree,50)
# 	sparse.diff.otu.2.plot.data[i,] <- unlist(data.frame(t(sparse.diff.otu.2.reps[[i]]$effect)))
# 	sparse.diff.2.dist[i,] <- c(sparse.diff.2.reps[[i]]$meanDist.SD[[1]]$meanDist,sparse.diff.2.reps[[i]]$meanDist.SD[[2]]$meanDist,sparse.diff.2.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	sparse.diff.2.sd[i,] <- c(sparse.diff.2.reps[[i]]$meanDist.SD[[1]]$error,sparse.diff.2.reps[[i]]$meanDist.SD[[2]]$error,sparse.diff.2.reps[[i]]$meanDist.SD[[3]]$error)
# }

# #med/high
# sparse.diff.otu.3.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# sparse.diff.otu.3.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.diff.otu.3.plot.data) <- plotDataTextLabels

# sparse.diff.3.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.diff.3.dist) <- plotDataTextLabels
# sparse.diff.3.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(sparse.diff.3.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	sparse.diff.otu.3.reps[[i]] <- runMixedReplicate(sparse.diff.otu.0001,sparse.diff.otu.00001,high.groups,high.groups,high.tree,50)
# 	sparse.diff.otu.3.plot.data[i,] <- unlist(data.frame(t(sparse.diff.otu.3.reps[[i]]$effect)))
# 	sparse.diff.3.dist[i,] <- c(sparse.diff.3.reps[[i]]$meanDist.SD[[1]]$meanDist,sparse.diff.3.reps[[i]]$meanDist.SD[[2]]$meanDist,sparse.diff.3.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	sparse.diff.3.sd[i,] <- c(sparse.diff.3.reps[[i]]$meanDist.SD[[1]]$error,sparse.diff.3.reps[[i]]$meanDist.SD[[2]]$error,sparse.diff.3.reps[[i]]$meanDist.SD[[3]]$error)
# }

# pdf("sparsityDifferenceTestPlots.pdf")
# par(mar=c(13, 4, 4, 2) + 0.1)
# par(cex.lab=1.3)
# par(cex.main=1.5)
# stripchart(sparse.diff.otu.1.plot.data,vertical=TRUE,main="Sparsity filter at 0.1% vs 0.01%",group.names=colnames(sparse.diff.otu.1.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# stripchart(sparse.diff.otu.2.plot.data,vertical=TRUE,main="Sparsity filter at 0.1% vs 0.001%",group.names=colnames(sparse.diff.otu.2.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# stripchart(sparse.diff.otu.3.plot.data,vertical=TRUE,main="Sparsity filter at 0.01% vs 0.001%",group.names=colnames(sparse.diff.otu.3.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# par(originalPar)

# par(mar=c(13, 6, 4, 2) + 0.1)

# meanDist <- unlist(lapply(sparse.diff.1.dist,mean))
# meanSd <- unlist(lapply(sparse.diff.1.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sparsity filter at 0.1% vs 0.01%")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# meanDist <- unlist(lapply(sparse.diff.2.dist,mean))
# meanSd <- unlist(lapply(sparse.diff.2.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sparsity filter at 0.1% vs 0.001%")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# meanDist <- unlist(lapply(sparse.diff.3.dist,mean))
# meanSd <- unlist(lapply(sparse.diff.3.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sparsity filter at 0.01% vs 0.001%")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# par(originalPar)

# dev.off()


# SHANNON DIVERSITY TEST

#note that average diversity isn't the same in the different conditions
#	medians are 6.216 for saliva, 5.624 for stool
high.otu.diversity <- diversity(high.otu)

#low/high diversity cutoffs set so that there are at least 10 samples in each condition
#	less samples -> extra replicates
high.diversity.indices <- which(high.otu.diversity > 6)
low.diversity.indices <- which(high.otu.diversity < 5.7)
high.diversity <- high.otu[high.diversity.indices,]
low.diversity <- high.otu[low.diversity.indices,]
high.diversity.groups <- high.groups[high.diversity.indices]
low.diversity.groups <- high.groups[low.diversity.indices]


diversity <- list()

diversity$low <- low.diversity
diversity$high <- high.diversity

diversity.groups <- list()
diversity.groups$low <- low.diversity.groups
diversity.groups$high <- high.diversity.groups

diversity.tree <- list()
diversity.tree[[1]] <- diversity.tree[[2]] <- high.tree

comparisonTitleList <- c("Shannon diversity < 5.7","Shannon diversity > 6")
fileName <- "diversityTestPlots.pdf"
nSamples <- 10

compare(runReplicate, diversity, diversity.groups, selfComparisonList2, diversity.tree, comparisonTitleList, fileName, nSamples)


# #low diversity
# low.diversity.otu.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# low.diversity.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.diversity.plot.data) <- plotDataTextLabels

# low.diversity.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.diversity.dist) <- plotDataTextLabels
# low.diversity.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.diversity.sd) <- plotDataTextLabels

# for (i in 1:extraReplicates) {
# 	low.diversity.otu.reps[[i]] <- runReplicate(low.diversity,low.diversity.groups,high.tree,10)
# 	low.diversity.plot.data[i,] <- unlist(data.frame(t(low.diversity.otu.reps[[i]]$effect)))
# 	low.diversity.dist[i,] <- c(low.diversity.reps[[i]]$meanDist.SD[[1]]$meanDist,low.diversity.reps[[i]]$meanDist.SD[[2]]$meanDist,low.diversity.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	low.diversity.sd[i,] <- c(low.diversity.reps[[i]]$meanDist.SD[[1]]$error,low.diversity.reps[[i]]$meanDist.SD[[2]]$error,low.diversity.reps[[i]]$meanDist.SD[[3]]$error)
# }

# #high diversity
# high.diversity.otu.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# high.diversity.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(high.diversity.plot.data) <- plotDataTextLabels

# high.diversity.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(high.diversity.dist) <- plotDataTextLabels
# high.diversity.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(high.diversity.sd) <- plotDataTextLabels

# for (i in 1:extraReplicates) {
# 	high.diversity.otu.reps[[i]] <- runReplicate(high.diversity,high.diversity.groups,high.tree,10)
# 	high.diversity.plot.data[i,] <- unlist(data.frame(t(high.diversity.otu.reps[[i]]$effect)))
# 	high.diversity.dist[i,] <- c(high.diversity.reps[[i]]$meanDist.SD[[1]]$meanDist,high.diversity.reps[[i]]$meanDist.SD[[2]]$meanDist,high.diversity.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	high.diversity.sd[i,] <- c(high.diversity.reps[[i]]$meanDist.SD[[1]]$error,high.diversity.reps[[i]]$meanDist.SD[[2]]$error,high.diversity.reps[[i]]$meanDist.SD[[3]]$error)
# }

# pdf("diversityTestPlots.pdf")
# par(mar=c(13, 4, 4, 2) + 0.1)
# par(cex.lab=1.3)
# par(cex.main=1.5)
# stripchart(low.diversity.plot.data,vertical=TRUE,main="Shannon diversity < 5.7",group.names=colnames(low.diversity.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# stripchart(high.diversity.plot.data,vertical=TRUE,main="Shannon diversity > 6",group.names=colnames(high.diversity.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# par(originalPar)

# par(mar=c(13, 6, 4, 2) + 0.1)

# meanDist <- unlist(lapply(low.diversity.dist,mean))
# meanSd <- unlist(lapply(low.diversity.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Shannon diversity < 5.7")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# meanDist <- unlist(lapply(high.diversity.dist,mean))
# meanSd <- unlist(lapply(high.diversity.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Shannon diversity > 6")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# par(originalPar)

# dev.off()


# SHANNON DIVERSITY DIFFERENCE TEST

comparisonTitleList <- c("Saliva diversity < 5.7 vs. stool diversity > 6.1","Stool diversity < 5.7 vs. saliva diversity > 6.1")
fileName <- "diversityDifferenceTestPlots.pdf"
nSamples <- 10

compare(runMixedReplicate, diversity, diversity.groups, mixedComparisonList2, diversity.tree, comparisonTitleList, fileName, nSamples)



# #low/high diversity (saliva vs. stool)
# low.high.diversity.diff.otu.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# low.high.diversity.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.high.diversity.diff.plot.data) <- plotDataTextLabels

# low.high.diversity.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.high.diversity.dist) <- plotDataTextLabels
# low.high.diversity.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(low.high.diversity.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	low.high.diversity.diff.otu.reps[[i]] <- runMixedReplicate(low.diversity,high.diversity,low.diversity.groups,high.diversity.groups,high.tree,10)
# 	low.high.diversity.diff.plot.data[i,] <- unlist(data.frame(t(low.high.diversity.diff.otu.reps[[i]]$effect)))
# 	low.high.diversity.dist[i,] <- c(low.high.diversity.reps[[i]]$meanDist.SD[[1]]$meanDist,low.high.diversity.reps[[i]]$meanDist.SD[[2]]$meanDist,low.high.diversity.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	low.high.diversity.sd[i,] <- c(low.high.diversity.reps[[i]]$meanDist.SD[[1]]$error,low.high.diversity.reps[[i]]$meanDist.SD[[2]]$error,low.high.diversity.reps[[i]]$meanDist.SD[[3]]$error)
# }

# #high/low diversity (saliva vs. stool)
# high.low.diversity.diff.otu.reps <- list()
# #columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
# high.low.diversity.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
# colnames(high.low.diversity.diff.plot.data) <- plotDataTextLabels

# high.low.diversity.dist <- data.frame(matrix(nrow=5,ncol=9))
# colnames(high.low.diversity.dist) <- plotDataTextLabels
# high.low.diversity.sd <- data.frame(matrix(nrow=5,ncol=9))
# colnames(high.low.diversity.sd) <- plotDataTextLabels

# for (i in 1:replicates) {
# 	high.low.diversity.diff.otu.reps[[i]] <- runMixedReplicate(high.diversity,low.diversity,high.diversity.groups,low.diversity.groups,high.tree,10)
# 	high.low.diversity.diff.plot.data[i,] <- unlist(data.frame(t(high.low.diversity.diff.otu.reps[[i]]$effect)))
# 	high.low.diversity.dist[i,] <- c(high.low.diversity.reps[[i]]$meanDist.SD[[1]]$meanDist,high.low.diversity.reps[[i]]$meanDist.SD[[2]]$meanDist,high.low.diversity.reps[[i]]$meanDist.SD[[3]]$meanDist)
# 	high.low.diversity.sd[i,] <- c(high.low.diversity.reps[[i]]$meanDist.SD[[1]]$error,high.low.diversity.reps[[i]]$meanDist.SD[[2]]$error,high.low.diversity.reps[[i]]$meanDist.SD[[3]]$error)
# }

# pdf("diversityDifferenceTestPlots.pdf")
# par(mar=c(13, 4, 4, 2) + 0.1)
# par(cex.lab=1.3)
# par(cex.main=1.5)
# stripchart(low.high.diversity.diff.plot.data,vertical=TRUE,main="Saliva diversity < 5.7 vs. stool diversity > 6.1",group.names=colnames(low.high.diversity.diff.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# stripchart(high.low.diversity.diff.plot.data,vertical=TRUE,main="Stool diversity < 5.7 vs. saliva diversity > 6.1",group.names=colnames(high.low.diversity.diff.plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
# par(originalPar)

# par(mar=c(13, 6, 4, 2) + 0.1)

# meanDist <- unlist(lapply(low.high.diversity.dist,mean))
# meanSd <- unlist(lapply(low.high.diversity.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Saliva diversity < 5.7 vs. tool diversity > 6.1")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# meanDist <- unlist(lapply(high.low.diversity.dist,mean))
# meanSd <- unlist(lapply(high.low.diversity.sd,mean))
# myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Stool diversity < 5.7 vs. saliva diversity > 6.1")
# segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
# segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

# par(originalPar)

# dev.off()

# par(originalPar)








