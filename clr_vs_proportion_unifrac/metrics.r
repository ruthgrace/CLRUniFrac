#make overlap matrix

library(vegan)

getOverlap <- function (otu.tab) {	

	# Convert into CLR
	otu.tab <- as.matrix(otu.tab)
	otu.tab.original <- otu.tab
	otu.tab[otu.tab==0] <- 0.5
	otu.tab.clr <- apply(otu.tab, 1, function(x){log2(x) - mean(log2(x))})
	numCols <- ncol(otu.tab.clr)

	overlap <- matrix(nrow=numCols,ncol=numCols)

	minval <- min(otu.tab.clr)

	otu.tab.clr.positive <- otu.tab.clr - minval

	for(i in 1:numCols) {
		for (j in i:numCols) {
			compare <- data.frame(otu.tab.clr.positive[,i],otu.tab.clr.positive[,j])
			overlap[i,j] <- sum(apply(compare,1,min)) / sum(apply(compare,1,max))
			overlap[j,i] <- overlap[i,j]
		}
	}

	return(overlap)
}


averageReadCount <- function(otu.tab) {
	otu.tab <- as.matrix(otu.tab)
	totalReadCount <- rowSums(otu.tab,na.rm=TRUE)
	numSamples <- nrow(otu.tab)

	avgReadCount <- matrix(nrow=numSamples,ncol=numSamples)

	return(averageGeneric(numSamples,totalReadCount))
}

averageGeneric <- function(numSamples,totalReadCount) {
	avgReadCount <- matrix(nrow=numSamples,ncol=numSamples)

	for(i in 1:numSamples) {
		for (j in i:numSamples) {
			avgReadCount[i,j] <- mean(c(totalReadCount[i],totalReadCount[j]))
			avgReadCount[j,i] <- avgReadCount[i,j]
		}
	}
	return(avgReadCount)
}

maxGeneric <- function(numSamples,totalReadCount) {
	avgReadCount <- matrix(nrow=numSamples,ncol=numSamples)

	for(i in 1:numSamples) {
		for (j in i:numSamples) {
			avgReadCount[i,j] <- max(c(totalReadCount[i],totalReadCount[j]))
			avgReadCount[j,i] <- avgReadCount[i,j]
		}
	}
	return(avgReadCount)
}

minGeneric <- function(numSamples,totalReadCount) {
	avgReadCount <- matrix(nrow=numSamples,ncol=numSamples)

	for(i in 1:numSamples) {
		for (j in i:numSamples) {
			avgReadCount[i,j] <- min(c(totalReadCount[i],totalReadCount[j]))
			avgReadCount[j,i] <- avgReadCount[i,j]
		}
	}
	return(avgReadCount)
}

printSeparation <- function(ruthClrUnifrac.pcoa,gUnifrac.pcoa,eUnifrac.pcoa,condition1,condition2,groups) {

	print(paste("COMPARING",condition1,"and",condition2))

	#calculate separation clrunifrac
	group1.indices <- which(groups==condition1)
	group2.indices <- which(groups==condition2)

	sd.1 <- max(sd(ruthClrUnifrac.pcoa$vectors[group1.indices,1]),sd(ruthClrUnifrac.pcoa$vectors[group2.indices,1]))
	sd.2 <- max(sd(ruthClrUnifrac.pcoa$vectors[group1.indices,2]),sd(ruthClrUnifrac.pcoa$vectors[group2.indices,2]))
	sd.3 <- max(sd(ruthClrUnifrac.pcoa$vectors[group1.indices,3]),sd(ruthClrUnifrac.pcoa$vectors[group2.indices,3]))

	group1.1 <- mean(ruthClrUnifrac.pcoa$vectors[group1.indices,1])/sd.1
	group2.1 <- mean(ruthClrUnifrac.pcoa$vectors[group2.indices,1])/sd.1

	group1.2 <- mean(ruthClrUnifrac.pcoa$vectors[group1.indices,2])/sd.2
	group2.2 <- mean(ruthClrUnifrac.pcoa$vectors[group2.indices,2])/sd.2

	group1.3 <- mean(ruthClrUnifrac.pcoa$vectors[group1.indices,3])/sd.3
	group2.3 <- mean(ruthClrUnifrac.pcoa$vectors[group2.indices,3])/sd.3

	group1.12 <- sqrt(group1.1^2 + group1.2^2)
	group2.12 <- sqrt(group2.1^2 + group2.2^2)

	group1.123 <- sqrt(group1.12^2 + group1.3^2)
	group2.123 <- sqrt(group2.12^2 + group2.3^2)



	print(paste("CLRUniFrac 1st component separation: ",abs(group1.1-group2.1)))
	print(paste("CLRUniFrac 1st and 2nd component separation: ",abs(group1.12-group2.12)))
	print(paste("CLRUniFrac 1,2,3 component separation: ",abs(group1.123-group2.123)))


	#calculate separation gunifrac
	group1.1 <- mean(gUnifrac.pcoa$vectors[group1.indices,1])
	group2.1 <- mean(gUnifrac.pcoa$vectors[group2.indices,1])

	group1.2 <- mean(gUnifrac.pcoa$vectors[group1.indices,2])
	group2.2 <- mean(gUnifrac.pcoa$vectors[group2.indices,2])

	group1.3 <- mean(gUnifrac.pcoa$vectors[group1.indices,3])
	group2.3 <- mean(gUnifrac.pcoa$vectors[group2.indices,3])

	group1.12 <- sqrt(group1.1^2 + group1.2^2)
	group2.12 <- sqrt(group2.1^2 + group2.2^2)

	group1.123 <- sqrt(group1.12^2 + group1.3^2)
	group2.123 <- sqrt(group2.12^2 + group2.3^2)

	print(paste("gUniFrac 1st component separation: ",abs(group1.1-group2.1)))
	print(paste("gUniFrac 1st and 2nd component separation: ",abs(group1.12-group2.12)))
	print(paste("gUniFrac 1,2,3 component separation: ",abs(group1.123-group2.123)))


	#calculate separation eunifrac
	group1.1 <- mean(eUnifrac.pcoa$vectors[group1.indices,1])
	group2.1 <- mean(eUnifrac.pcoa$vectors[group2.indices,1])

	group1.2 <- mean(eUnifrac.pcoa$vectors[group1.indices,2])
	group2.2 <- mean(eUnifrac.pcoa$vectors[group2.indices,2])

	group1.3 <- mean(eUnifrac.pcoa$vectors[group1.indices,3])
	group2.3 <- mean(eUnifrac.pcoa$vectors[group2.indices,3])

	group1.12 <- sqrt(group1.1^2 + group1.2^2)
	group2.12 <- sqrt(group2.1^2 + group2.2^2)

	group1.123 <- sqrt(group1.12^2 + group1.3^2)
	group2.123 <- sqrt(group2.12^2 + group2.3^2)

	print(paste("eUniFrac 1st component separation: ",abs(group1.1-group2.1)))
	print(paste("eUniFrac 1st and 2nd component separation: ",abs(group1.12-group2.12)))
	print(paste("eUniFrac 1,2,3 component separation: ",abs(group1.123-group2.123)))

}

getShannonDiversityDiffMat <- function(otu) {
	diversityList <- diversity(otu)
	#put diversity into rows
	diversityList <- t(diversityList)
	diversityList <- t(diversityList)
	distMat <- dist(diversityList,method="manhattan",diag=TRUE,upper=TRUE)
	return(distMat)
}

getAvgShannonDiversity <- function(otu) {
	diversityList <- diversity(otu)
	return(averageGeneric(nrow(otu),diversityList))
}

getMinShannonDiversity <- function(otu) {
	diversityList <- diversity(otu)
	return(minGeneric(nrow(otu),diversityList))
}

getMaxShannonDiversity <- function(otu) {
	diversityList <- diversity(otu)
	return(maxGeneric(nrow(otu),diversityList))
}

getDataSetSep <- function(otu,groups,tree) {
	#take in otu table, rows are samples, cols are OTUs
	# return list of 1) unweighted UniFrac 2) weighted UniFrac 3) iUniFrac

	unifrac <- GUniFrac(otu, brazil.tree, alpha = c(1))
	uwUnifrac <- unifrac$unifrac[,,1]
	wUnifrac <- unifrac$unifrac[,,3]
	eUnifrac <- InformationUniFrac(otu, brazil.tree, alpha = c(1))$unifrac[,,1]

	uwUnifrac.pcoa <- pcoa(uwUnifrac)
	wUnifrac.pcoa <- pcoa(wUnifrac)
	eUnifrac.pcoa <- pcoa(eUnifrac)

	uwUnifrac.sep <- getPCoASep(uwUnifrac.pcoa)
	wUnifrac.sep <- getPCoASep(wUnifrac.pcoa)
	eUnifrac.sep <- getPCoASep(eUnifrac.pcoa)

	returnList <- data.frame(uwUnifrac.sep,wUnifrac.sep,eUnifrac.sep)
	colnames(returnList) <- c("uwUnifrac","wUnifrac","eUnifrac")
	return(returnList)
}

getPCoASep <- function(pcoa,groups) {
	#given pcoa and metadata
	# returns separations on axis 1, 1&2, 1&2&3
	group1.1 <- pcoa$vector[which(groups==levels(groups)[1]),1]
	group2.1 <- pcoa$vector[which(groups==levels(groups)[2]),1]
	diff.1 <- abs(mean(group1.1) - mean(group2.1))/sd(pcoa$vector[,1])

	group1.2 <- pcoa$vector[which(groups==levels(groups)[1]),2]
	group2.2 <- pcoa$vector[which(groups==levels(groups)[2]),2]
	diff.2 <- abs(mean(group1.1) - mean(group2.1))/sd(pcoa$vector[,2])

	group1.3 <- pcoa$vector[which(groups==levels(groups)[1]),3]
	group2.3 <- pcoa$vector[which(groups==levels(groups)[2]),3]
	diff.3 <- abs(mean(group1.1) - mean(group2.1))/sd(pcoa$vector[,3])

	diff.12 <- sqrt(diff.1^2 + diff.2^2)
	diff.123 <- sqrt(diff.12^2 + diff.3^2)

	returnList <- data.frame(c(diff.1,diff.12,diff.123))
	rownames(returnList) <- c("separationOn1,separationOn12,separationOn123")
	returnList <- t(returnList)
	return(returnList)
}