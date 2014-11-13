#make overlap matrix

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

	for(i in 1:numSamples) {
		for (j in i:numSamples) {
			avgReadCount[i,j] <- mean(c(totalReadCount[i],totalReadCount[j]))
			avgReadCount[j,i] <- avgReadCount[i,j]
		}
	}
	return(avgReadCount)
}
