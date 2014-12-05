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

averageGeneric <- function(numSamples,matrix) {
	avgReadCount <- matrix(nrow=numSamples,ncol=numSamples)

	for(i in 1:numSamples) {
		for (j in i:numSamples) {
			avgReadCount[i,j] <- mean(c(totalReadCount[i],totalReadCount[j]))
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
