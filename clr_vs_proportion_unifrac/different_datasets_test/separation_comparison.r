# consider testing with different rarification methods later ?

options(error=recover)

library(ape)
library(phangorn)
library(vegan)

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

runReplicate <- function(otu,groups,tree) {
	
	nSamples <- 50
	
	#sample 50 samples from condition 1
	group1.indices <- which(groups==levels(groups)[1])
	group1.rand <- otu[as.integer(sample(group1.indices,nSamples,replace=FALSE)),]
	
	#sample 50 samples from condition 2
	group2.indices <- which(groups==levels(groups)[2])
	group2.rand <- otu[as.integer(sample(group2.indices,nSamples,replace=FALSE)),]
	
	#concatenate
	data <- rbind(group1.rand,group2.rand)

	#make groups
	newGroups <- c(rep(levels(groups)[1],50),rep(levels(groups)[2],50))
	return(getDataSetSep(data,newGroups,tree))
}

runMixedReplicate <- function(otu1,otu2,groups1,groups2,tree) {
	
	nSamples <- 50
	
	#sample 50 samples from condition 1,group1
	group1.indices <- which(groups1==levels(groups1)[1])
	group1.rand <- otu1[as.integer(sample(group1.indices,nSamples,replace=FALSE)),]
	
	#sample 50 samples from condition 2,group2
	group2.indices <- which(groups2==levels(groups2)[2])
	group2.rand <- otu2[as.integer(sample(group2.indices,nSamples,replace=FALSE)),]
	
	#concatenate
	data <- rbind(group1.rand,group2.rand)

	#make groups
	newGroups <- c(rep(levels(groups)[1],50),rep(levels(groups)[2],50))
	return(getDataSetSep(data,newGroups,tree))
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
plotSeqDepthDataColNames <- paste(pcoaLabels, unifracLabels, sep = ".")

sparsityLabels <- c(rep("sparsity.001",3),rep("sparsity.0001",3),rep("sparsity.00001",3))
plotSparsityDataColNames <- paste(sparsityLabels, unifracLabels, sep = ".")

diversityLabels <- c(rep("low.diversity",3),rep("high.diversity",3))
shortUnifracLabels <- rep(c("uwUnifrac","wUnifrac","iUnifrac"),2)
plotDiversityDataColNames <- paste(diversityLabels, shortUnifracLabels, sep = ".")


# SEQUENCING DEPTH TEST
#low
low.seq.depth.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
low.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(low.seq.depth.plot.data) <- plotSeqDepthDataColNames
for (i in 1:replicates) {
	low.seq.depth.reps[[i]] <- runReplicate(low.otu,low.groups,low.tree)
	low.seq.depth.plot.data[i,] <- unlist(data.frame(t(low.seq.depth.reps[[i]])))
}
#med
med.seq.depth.reps <- list()
med.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(med.seq.depth.plot.data) <- plotSeqDepthDataColNames
for (i in 1:replicates) {
	med.seq.depth.reps[[i]] <- runReplicate(med.otu,med.groups,med.tree)
	med.seq.depth.plot.data[i,] <- unlist(data.frame(t(med.seq.depth.reps[[i]])))
}
#high
high.seq.depth.reps <- list()
high.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(high.seq.depth.plot.data) <- plotSeqDepthDataColNames
for (i in 1:replicates) {
	high.seq.depth.reps[[i]] <- runReplicate(high.otu,high.groups,high.tree)
	high.seq.depth.plot.data[i,] <- unlist(data.frame(t(high.seq.depth.reps[[i]])))
}

#plot
stripchart(data.frame(),pch=19,col=transparentdarkorchid)

# SEQUENCING DEPTH DIFFERENCE TEST
#low/med
low.seq.depth.diff.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
low.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(low.seq.depth.diff.plot.data) <- plotSeqDepthDataColNames
for (i in 1:replicates) {
	low.seq.depth.diff.reps[[i]] <- runMixedReplicate(low.otu,med.otu,low.groups,med.groups,low.tree)
	low.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(low.seq.depth.diff.reps[[i]])))
}
#med/high
med.seq.depth.diff.reps <- list()
med.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(med.seq.depth.diff.plot.data) <- plotSeqDepthDataColNames
for (i in 1:replicates) {
	med.seq.depth.diff.reps[[i]] <- runMixedReplicate(med.otu,high.otu,med.groups,high.groups,med.tree)
	med.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(med.seq.depth.diff.reps[[i]])))
}
#low/high
high.seq.depth.diff.reps <- list()
high.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(high.seq.depth.diff.plot.data) <- plotSeqDepthDataColNames
for (i in 1:replicates) {
	high.seq.depth.diff.reps[[i]] <- runMixedReplicate(low.otu,high.otu,low.groups,high.groups,high.tree)
	high.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(high.seq.depth.diff.reps[[i]])))
}

# SPARSITY TEST
#remove OTUs rarer than a thresh hold throughout all samples
high.otu.sum <- apply(high.otu,2,sum)
high.total.sum <- sum(high.otu)

#0.1% sparsity filter
sparse.otu.001 <- high.otu[,(which(high.otu.sum >= (0.001*high.total.sum)))]
sparse.otu.001.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.otu.001.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.otu.001.plot.data) <- plotSparsityDataColNames
for (i in 1:replicates) {
	sparse.otu.001.reps[[i]] <- runReplicate(sparse.otu.001,high.groups,high.tree)
	sparse.otu.001.plot.data[i,] <- unlist(data.frame(t(sparse.otu.001.reps[[i]])))
}

#0.01% sparsity filter
sparse.otu.0001 <- high.otu[,(which(high.otu.sum >= (0.0001*high.total.sum)))]
sparse.otu.0001.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.otu.0001.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.otu.0001.plot.data) <- plotSparsityDataColNames
for (i in 1:replicates) {
	sparse.otu.0001.reps[[i]] <- runReplicate(sparse.otu.0001,high.groups,high.tree)
	sparse.otu.0001.plot.data[i,] <- unlist(data.frame(t(sparse.otu.0001.reps[[i]])))
}

#0.001% sparsity filter
sparse.otu.00001 <- high.otu[,(which(high.otu.sum >= (0.00001*high.total.sum)))]
sparse.otu.00001.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.otu.00001.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.otu.00001.plot.data) <- plotSparsityDataColNames
for (i in 1:replicates) {
	sparse.otu.00001.reps[[i]] <- runReplicate(sparse.otu.00001,high.groups,high.tree)
	sparse.otu.00001.plot.data[i,] <- unlist(data.frame(t(sparse.otu.00001.reps[[i]])))
}

# SPARSITY DIFFERENCE TEST

#low/med
sparse.diff.otu.1.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.diff.otu.1.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.diff.otu.1.plot.data) <- plotSparsityDataColNames
for (i in 1:replicates) {
	sparse.diff.otu.1.reps[[i]] <- runMixedReplicate(sparse.diff.otu.001,sparse.diff.otu.0001,high.groups,high.groups,high.tree)
	sparse.diff.otu.1.plot.data[i,] <- unlist(data.frame(t(sparse.diff.otu.1.reps[[i]])))
}

#low/high
sparse.diff.otu.2.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.diff.otu.2.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.diff.otu.2.plot.data) <- plotSparsityDataColNames
for (i in 1:replicates) {
	sparse.diff.otu.2.reps[[i]] <- runMixedReplicate(sparse.diff.otu.001,sparse.diff.otu.00001,high.groups,high.groups,high.tree)
	sparse.diff.otu.2.plot.data[i,] <- unlist(data.frame(t(sparse.diff.otu.2.reps[[i]])))
}

#med/high
sparse.diff.otu.3.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
sparse.diff.otu.3.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(sparse.diff.otu.3.plot.data) <- plotSparsityDataColNames
for (i in 1:replicates) {
	sparse.diff.otu.3.reps[[i]] <- runMixedReplicate(sparse.diff.otu.0001,sparse.diff.otu.00001,high.groups,high.groups,high.tree)
	sparse.diff.otu.3.plot.data[i,] <- unlist(data.frame(t(sparse.diff.otu.3.reps[[i]])))
}


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
colnames(low.diversity.plot.data) <- plotDiversityDataColNames
for (i in 1:replicates) {
	low.diversity.otu.reps[[i]] <- runReplicate(low.diversity,low.diversity.groups,high.tree)
	low.diversity.plot.data[i,] <- unlist(data.frame(t(low.diversity.otu.reps[[i]])))
}

#high diversity
high.diversity.otu.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
high.diversity.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(high.diversity.plot.data) <- plotDiversityDataColNames
for (i in 1:replicates) {
	high.diversity.otu.reps[[i]] <- runReplicate(high.diversity,high.diversity.groups,high.tree)
	high.diversity.plot.data[i,] <- unlist(data.frame(t(high.diversity.otu.reps[[i]])))
}


# SHANNON DIVERSITY DIFFERENCE TEST

#low/high diversity (saliva vs. stool)

#high/low diversity (saliva vs. stool)













