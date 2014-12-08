# consider testing with different rarification methods later ?

options(error=recover)

library(ape)
library(phangorn)

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

darkorchid <- col2rgb("darkorchid4")
transparentdarkorchid <- rgb(darkorchid[1]/255,darkorchid[2]/255,darkorchid[3]/255,0.3)

pcoaLabels <- c(rep("pcoa1",3),rep("pcoa12",3),rep("pcoa123",3))
unifracLabels <- rep(c("uwUnifrac","wUnifrac","iUnifrac"),3)
plotDataColNames <- paste(pcoaLabels, unifracLabels, sep = ".")

# SEQUENCING DEPTH TEST
#low
low.seq.depth.reps <- list()
#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
low.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(low.seq.depth.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	low.seq.depth.reps[[i]] <- runReplicate(low.otu,low.groups,low.tree)
	low.seq.depth.plot.data[i,] <- unlist(data.frame(t(low.seq.depth.reps[[i]])))
}
#med
med.seq.depth.reps <- list()
med.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(med.seq.depth.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	med.seq.depth.reps[[i]] <- runReplicate(med.otu,med.groups,med.tree)
	med.seq.depth.plot.data[i,] <- unlist(data.frame(t(med.seq.depth.reps[[i]])))
}
#high
high.seq.depth.reps <- list()
high.seq.depth.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(high.seq.depth.plot.data) <- plotDataColNames
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
colnames(low.seq.depth.diff.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	low.seq.depth.diff.reps[[i]] <- runMixedReplicate(low.otu,med.otu,low.groups,med.groups,low.tree)
	low.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(low.seq.depth.diff.reps[[i]])))
}
#med/high
med.seq.depth.diff.reps <- list()
med.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(med.seq.depth.diff.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	med.seq.depth.diff.reps[[i]] <- runMixedReplicate(med.otu,high.otu,med.groups,high.groups,med.tree)
	med.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(med.seq.depth.diff.reps[[i]])))
}
#low/high
high.seq.depth.diff.reps <- list()
high.seq.depth.diff.plot.data <- data.frame(matrix(nrow=5,ncol=9))
colnames(high.seq.depth.diff.plot.data) <- plotDataColNames
for (i in 1:replicates) {
	high.seq.depth.diff.reps[[i]] <- runMixedReplicate(low.otu,high.otu,low.groups,high.groups,high.tree)
	high.seq.depth.diff.plot.data[i,] <- unlist(data.frame(t(high.seq.depth.diff.reps[[i]])))
}



# SPARSITY TEST

# SPARSITY DIFFERENCE TEST

#sparsity filter
#remove all OTUs for which the minimum count is < 30
# mouth.otu.min <- apply(mouth.otu,2,min)
# mouth.otu <- mouth.otu[,which(mouth.otu.min) >= 30)]
# mouth.original <- mouth.original[,which(mouth.otu.min) >= 30)]


# SHANNON DIVERSITY TEST

# SHANNON DIVERSITY DIFFERENCE TEST















