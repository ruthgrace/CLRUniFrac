# consider testing with different rarification methods later ?

options(error=recover)

library(ape)
library(phangorn)
library(vegan)

originalPar <- par()
par(xaxt="n")

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

source("../metrics.r")

replicates <- 5
extraReplicates <- 10

darkorchid <- col2rgb("darkorchid4")
transparentdarkorchid <- rgb(darkorchid[1]/255,darkorchid[2]/255,darkorchid[3]/255,0.3)

pcoaLabels <- c(rep("pcoa1",3),rep("pcoa12",3),rep("pcoa123",3))
unifracLabels <- rep(c("uwUnifrac","wUnifrac","iUnifrac"),3)
plotDataColNames <- paste(pcoaLabels, unifracLabels, sep = ".")

# sparsityLabels <- c(rep("sparsity.001",3),rep("sparsity.0001",3),rep("sparsity.00001",3))
# plotSparsityDataColNames <- paste(sparsityLabels, unifracLabels, sep = ".")

# diversityLabels <- c(rep("low.diversity",3),rep("high.diversity",3))
# shortUnifracLabels <- rep(c("uwUnifrac","wUnifrac","iUnifrac"),2)
# plotDiversityDataColNames <- paste(diversityLabels, shortUnifracLabels, sep = ".")


# SEQUENCING DEPTH TEST
#low
low.seq.depth.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
low.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(low.seq.depth.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	low.seq.depth.reps[[i]] <- runReplicate(low.otu,low.groups,low.tree,50)
	low.seq.depth.plot.data[i,] <- unlist(data.frame(t(low.seq.depth.reps[[i]]$effect)))
}
#med
med.seq.depth.reps <- list()
med.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(med.seq.depth.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	med.seq.depth.reps[[i]] <- runReplicate(med.otu,med.groups,med.tree,50)
	med.seq.depth.plot.data[i,] <- unlist(data.frame(t(med.seq.depth.reps[[i]]$effect)))
}
#high
high.seq.depth.reps <- list()
high.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(high.seq.depth.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	high.seq.depth.reps[[i]] <- runReplicate(high.otu,high.groups,high.tree,50)
	high.seq.depth.plot.data[i,] <- unlist(data.frame(t(high.seq.depth.reps[[i]]$effect)))
}

#plot
pdf("sequencingDepthPlots.pdf")
stripchart(low.seq.depth.plot.data,vertical=TRUE,main="sequencing depth < 3000 reads/sample",group.names=colnames(low.seq.depth.plot.data),pch=19,col=transparentdarkorchid)
stripchart(med.seq.depth.plot.data,vertical=TRUE,main="sequencing depth 3000-6000 reads/sample",group.names=colnames(med.seq.depth.plot.data),pch=19,col=transparentdarkorchid)
stripchart(high.seq.depth.plot.data,vertical=TRUE,main="sequencing depth > 6000 reads/sample",group.names=colnames(high.seq.depth.plot.data),pch=19,col=transparentdarkorchid)
dev.off()



# SEQUENCING DEPTH DIFFERENCE TEST

#low/med
low.seq.depth.diff.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
low.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(low.seq.depth.diff.plot.data) <- plotDataColNames
addedOTUs <- addDisimilarOTUs(low.otu,med.otu,high.tree)
low.with.med.otu <- addedOTUs[[1]]
med.with.low.otu <- addedOTUs[[2]]
for (i in 1:replicates) {
	low.seq.depth.diff.reps[[i]] <- runMixedReplicate(low.with.med.otu,med.with.low.otu,low.groups,med.groups,high.tree,50)
	low.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(low.seq.depth.diff.reps[[i]]$effect)))
}
#med/high
med.seq.depth.diff.reps <- list()
med.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(med.seq.depth.diff.plot.data) <- plotDataColNames
addedOTUs <- addDisimilarOTUs(med.otu,high.otu,high.tree)
med.with.high.otu <- addedOTUs[[1]]
high.with.med.otu <- addedOTUs[[2]]
for (i in 1:replicates) {
	med.seq.depth.diff.reps[[i]] <- runMixedReplicate(med.with.high.otu,high.with.med.otu,med.groups,high.groups,high.tree,50)
	med.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(med.seq.depth.diff.reps[[i]]$effect)))
}
#low/high
high.seq.depth.diff.reps <- list()
high.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(high.seq.depth.diff.plot.data) <- plotDataColNames
addedOTUs <- addDisimilarOTUs(low.otu,high.otu,high.tree)
low.with.high.otu <- addedOTUs[[1]]
high.with.low.otu <- addedOTUs[[2]]
for (i in 1:replicates) {
	high.seq.depth.diff.reps[[i]] <- runMixedReplicate(low.with.high.otu,high.with.low.otu,low.groups,high.groups,high.tree,50)
	high.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(high.seq.depth.diff.reps[[i]]$effect)))
}

pdf("sequencingDepthDifferencePlots.pdf")
stripchart(low.seq.depth.diff.plot.data,vertical=TRUE,main="sequencing depth < 3000 reads/sample",group.names=colnames(low.seq.depth.diff.plot.data),pch=19,col=transparentdarkorchid)
stripchart(med.seq.depth.diff.plot.data,vertical=TRUE,main="sequencing depth 3000-6000 reads/sample",group.names=colnames(med.seq.depth.diff.plot.data),pch=19,col=transparentdarkorchid)
stripchart(high.seq.depth.diff.plot.data,vertical=TRUE,main="sequencing depth > 6000 reads/sample",group.names=colnames(high.seq.depth.diff.plot.data),pch=19,col=transparentdarkorchid)
dev.off()

# SPARSITY TEST
#remove OTUs rarer than a thresh hold throughout all samples
high.otu.sum <- apply(high.otu,2,sum)
high.total.sum <- sum(high.otu)

#0.1% sparsity filter
sparse.otu.001 <- high.otu[,(which(high.otu.sum >= (0.001*high.total.sum)))]
sparse.otu.001.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.otu.001.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.otu.001.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	sparse.otu.001.reps[[i]] <- runReplicate(sparse.otu.001,high.groups,high.tree,50)
	sparse.otu.001.plot.data[i,] <- unlist(data.frame(t(sparse.otu.001.reps[[i]]$effect)))
}

#0.01% sparsity filter
sparse.otu.0001 <- high.otu[,(which(high.otu.sum >= (0.0001*high.total.sum)))]
sparse.otu.0001.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.otu.0001.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.otu.0001.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	sparse.otu.0001.reps[[i]] <- runReplicate(sparse.otu.0001,high.groups,high.tree,50)
	sparse.otu.0001.plot.data[i,] <- unlist(data.frame(t(sparse.otu.0001.reps[[i]]$effect)))
}

#0.001% sparsity filter
sparse.otu.00001 <- high.otu[,(which(high.otu.sum >= (0.00001*high.total.sum)))]
sparse.otu.00001.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.otu.00001.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.otu.00001.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	sparse.otu.00001.reps[[i]] <- runReplicate(sparse.otu.00001,high.groups,high.tree,50)
	sparse.otu.00001.plot.data[i,] <- unlist(data.frame(t(sparse.otu.00001.reps[[i]]$effect)))
}



##### NEED TO MAKE AXIS LOOK RIGHT + plot other stuff

pdf("sparsityTestPlots.pdf")
stripchart(sparse.otu.001.plot.data,vertical=TRUE,main="sparsity filter at 0.1%",group.names=colnames(sparse.otu.001.plot.data),pch=19,col=transparentdarkorchid)
stripchart(sparse.otu.0001.plot.data,vertical=TRUE,main="sparsity filter at 0.01%",group.names=colnames(sparse.otu.0001.plot.data),pch=19,col=transparentdarkorchid)
stripchart(sparse.otu.00001.plot.data,vertical=TRUE,main="sparsity filter at 0.001%",group.names=colnames(sparse.otu.00001.plot.data),pch=19,col=transparentdarkorchid)
dev.off()

# SPARSITY DIFFERENCE TEST
sparse.diff.otu.001 <- high.otu
sparse.diff.otu.001[,(which(high.otu.sum < (0.001*high.total.sum)))] <- 0
sparse.diff.otu.0001 <- high.otu
sparse.diff.otu.0001[,(which(high.otu.sum < (0.0001*high.total.sum)))] <- 0
sparse.diff.otu.00001 <- high.otu
sparse.diff.otu.00001[,(which(high.otu.sum < (0.00001*high.total.sum)))] <- 0

#low/med
sparse.diff.otu.1.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.diff.otu.1.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.diff.otu.1.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	sparse.diff.otu.1.reps[[i]] <- runMixedReplicate(sparse.diff.otu.001,sparse.diff.otu.0001,high.groups,high.groups,high.tree,50)
	sparse.diff.otu.1.plot.data[i,] <- unlist(data.frame(t(sparse.diff.otu.1.reps[[i]]$effect)))
}

#low/high
sparse.diff.otu.2.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.diff.otu.2.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.diff.otu.2.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	sparse.diff.otu.2.reps[[i]] <- runMixedReplicate(sparse.diff.otu.001,sparse.diff.otu.00001,high.groups,high.groups,high.tree,50)
	sparse.diff.otu.2.plot.data[i,] <- unlist(data.frame(t(sparse.diff.otu.2.reps[[i]]$effect)))
}

#med/high
sparse.diff.otu.3.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.diff.otu.3.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.diff.otu.3.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	sparse.diff.otu.3.reps[[i]] <- runMixedReplicate(sparse.diff.otu.0001,sparse.diff.otu.00001,high.groups,high.groups,high.tree,50)
	sparse.diff.otu.3.plot.data[i,] <- unlist(data.frame(t(sparse.diff.otu.3.reps[[i]]$effect)))
}

pdf("sparsityDifferenceTestPlots.pdf")
stripchart(sparse.diff.otu.1.plot.data,vertical=TRUE,main="sparsity filter at 0.1% vs 0.01%",group.names=colnames(sparse.diff.otu.1.plot.data),pch=19,col=transparentdarkorchid)
stripchart(sparse.diff.otu.2.plot.data,vertical=TRUE,main="sparsity filter at 0.1% vs 0.001%",group.names=colnames(sparse.diff.otu.2.plot.data),pch=19,col=transparentdarkorchid)
stripchart(sparse.diff.otu.3.plot.data,vertical=TRUE,main="sparsity filter at 0.01% vs 0.001%",group.names=colnames(sparse.diff.otu.3.plot.data),pch=19,col=transparentdarkorchid)
dev.off()


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

#low diversity
low.diversity.otu.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
low.diversity.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(low.diversity.plot.data) <- plotDataColNames
for (i in 1:extraReplicates) {
	low.diversity.otu.reps[[i]] <- runReplicate(low.diversity,low.diversity.groups,high.tree,10)
	low.diversity.plot.data[i,] <- unlist(data.frame(t(low.diversity.otu.reps[[i]]$effect)))
}

#high diversity
high.diversity.otu.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
high.diversity.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(high.diversity.plot.data) <- plotDataColNames
for (i in 1:extraReplicates) {
	high.diversity.otu.reps[[i]] <- runReplicate(high.diversity,high.diversity.groups,high.tree,10)
	high.diversity.plot.data[i,] <- unlist(data.frame(t(high.diversity.otu.reps[[i]]$effect)))
}

pdf("diversityTestPlots.pdf")
stripchart(low.diversity.plot.data,vertical=TRUE,main="diversity < 5.7",group.names=colnames(low.diversity.plot.data),pch=19,col=transparentdarkorchid)
stripchart(high.diversity.plot.data,vertical=TRUE,main="diversity > 6",group.names=colnames(high.diversity.plot.data),pch=19,col=transparentdarkorchid)
dev.off()


# SHANNON DIVERSITY DIFFERENCE TEST

#low/high diversity (saliva vs. stool)
low.high.diversity.diff.otu.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
low.high.diversity.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(low.high.diversity.diff.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	low.high.diversity.diff.otu.reps[[i]] <- runMixedReplicate(low.diversity,high.diversity,low.diversity.groups,high.diversity.groups,high.tree,10)
	low.high.diversity.diff.plot.data[i,] <- unlist(data.frame(t(low.high.diversity.diff.otu.reps[[i]]$effect)))
}

#high/low diversity (saliva vs. stool)
high.low.diversity.diff.otu.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
high.low.diversity.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(high.low.diversity.diff.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	high.low.diversity.diff.otu.reps[[i]] <- runMixedReplicate(high.diversity,low.diversity,high.diversity.groups,low.diversity.groups,high.tree,10)
	high.low.diversity.diff.plot.data[i,] <- unlist(data.frame(t(high.low.diversity.diff.otu.reps[[i]]$effect)))
}

pdf("diversityDifferenceTestPlots.pdf")
stripchart(low.high.diversity.diff.plot.data,vertical=TRUE,main="saliva diversity < 5.7 vs. stool diversity > 6.1",group.names=colnames(low.high.diversity.diff.plot.data),pch=19,col=transparentdarkorchid)
stripchart(high.low.diversity.diff.plot.data,vertical=TRUE,main="stool diversity < 5.7 vs. saliva diversity > 6.1",group.names=colnames(high.low.diversity.diff.plot.data),pch=19,col=transparentdarkorchid)
dev.off()

par(originalPar)








