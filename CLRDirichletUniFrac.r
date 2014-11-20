##############################################################
# 
# Hacked by Ruth Grace Wong (ruthgracewong@gmail.com)
# to use weightings based on Centered Log Ratio transformed
# OTU count values (rather than proportional abundance)
#
# * doesn't do variance adjusted unifrac distances anymore.
#
#
#
# ORIGINAL PACKAGE INFORMATION (most of the code copied from here)
# 
# GUniFrac: Generalized UniFrac distances for comparing microbial
#						communities.
# Jun Chen (chenjun@mail.med.upenn.edu)
# Feb 24, 2012
#
# Reference: Jun Chen & Hongzhe (2012). Associating microbiome 
# composition with environmental covariates using generalized 
# UniFrac distances. (Submitted)
###############################################################

require(ade4)
require(ape)
require(vegan)



#dirichlet method - gets dirichlet distribution of reads per OTU, returns median
rdirichlet <- function (alpha)
#ISALIAS ddirichlet
#--------------------------------------------
{
	#128 replicates (median converges at about 64 replicates)
	n <- 128
	if(length(n) > 1) n <- length(n)
	if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
	n <- as.integer(n)
	if(n < 0) stop("integer(n) can not be negative in rtriang")

	if(is.vector(alpha)) alpha <- t(alpha)
	l <- dim(alpha)[2]
	x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  # Gere le recycling
	return(x / rowSums(x))
}


CLRDirichletUniFrac <- function (otu.tab, tree, alpha = c(0, 0.5, 1)) {
	# Calculate Generalized UniFrac distances. Unweighted and 
	# Variance-adjusted UniFrac distances will also be returned.
	#	
	# Args:
	#		otu.tab: OTU count table, row - n sample, column - q OTU
	#		tree: rooted phylogenetic tree of R class "phylo"
	#		alpha: parameter controlling weight on abundant lineages
	#
	# Returns:
	# 	unifracs: three dimensional array containing the generalized 
	#							UniFrac distances, unweighted UniFrac distance and 
	#							variance adjusted UniFrac distances. 
	#

	if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!")


	# Convert into CLR
	
	otu.tab <- data.matrix(otu.tab)
	#remove all OTUs that have no counts in all samples
	otu.tab <- otu.tab[which(rowSums(otu.tab)!=0),]
	# add prior of 0.5 -- zero counts are assumed to represent presence half the time
	otu.tab[otu.tab==0] <- 0.5
	otu.tab.mean <- apply(otu.tab, 1, function(x) {mean(log2(x))} )
	numOTUs <- rowSums(otu.tab!=0) # how many in each column is not zero

	#calculate proportions using the Dirichlet
	otu.tab.D <- t(apply(t(otu.tab), 2, function(x){rdirichlet(x)}))

	otu.tab.clr <- apply(otu.tab.D, 1, function(x){log2(x) - mean(log2(x))})
	otu.tab.clr <- t(otu.tab.clr)
	n <- nrow(otu.tab)

	# Construct the returning array
	if (is.null(rownames(otu.tab))) {
		rownames(otu.tab) <- paste("comm", 1:n, sep="_")
	}

	# d_UW: unweighted UniFrac, d_VAW: weighted UniFrac
	dimname3 <- c(paste("d", alpha, sep="_"), "d_UW", "d_VAW")
	unifracs <- array(NA, c(n, n, length(alpha) + 2),
				  dimnames=list(rownames(otu.tab), rownames(otu.tab), dimname3))
	for (i in 1:(length(alpha)+2)){
		for (j in 1:n){
			unifracs[j, j, i] <- 0
		}
	}	
	
	# Check OTU name consistency
	if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
		stop("The OTU table contains unknown OTUs! OTU names
					in the OTU table and the tree should match!" )
	}
	
	# Get the subtree if tree contains more OTUs
	absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
	if (length(absent) != 0) {
		tree <- drop.tip(tree, absent)
		warning("The tree has more OTU than the OTU table!")
	}
	
	# Reorder the otu.tab matrix if the OTU orders are different
	tip.label <- tree$tip.label
	otu.tab <- otu.tab[, tip.label]
	
	ntip <- length(tip.label)
	nbr <- nrow(tree$edge)	
	edge <- tree$edge
	edge2 <- edge[, 2]
	br.len <- tree$edge.length
	
	#  Accumulate OTU CLR measurements up the tree
	cum <- matrix(0, nbr, n)							# Branch abundance matrix
	effectiveNumSamples <- matrix(0, nbr, n)	#number of children
	effectiveNumSamples <- t(apply(effectiveNumSamples,1,function(x) numOTUs))
	clrcum <- matrix(0, nbr, n)							# clr cumulative
	geometricMean <- matrix(0, nbr, n)							# center of clr
	
	
	for (i in 1:ntip) {
		tip.loc <- which(edge2 == i)
		cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
		clrcum[tip.loc, ] <- otu.tab.clr[, i]
		geometricMean[tip.loc, ] <- otu.tab.mean
		node <- edge[tip.loc, 1]						# Assume the direction of edge 
		node.loc <- which(edge2 == node)
		while (length(node.loc)) {
			#calculate what the centered log ratio transform would be for this node if all the read counts under it were counted under one taxa
			oldMean <- geometricMean[node.loc, cum[node.loc, ]!=0 & otu.tab[, i]!=0]
			numSamples <- effectiveNumSamples[node.loc, cum[node.loc, ]!=0 & otu.tab[, i]!=0]
			x <- cum[node.loc, cum[node.loc, ]!=0 & otu.tab[, i]!=0]
			y <- otu.tab[cum[node.loc, ]!=0 & otu.tab[, i]!=0, i]
			newMean <- ( (oldMean*numSamples) - log2(x) - log2(y) + log2(x+y) )/numSamples-1
			geometricMean[node.loc, cum[node.loc, ]!=0 & otu.tab[, i]!=0] <- newMean
			geometricMean[node.loc,  geometricMean[node.loc,]==0]<- otu.tab.mean[geometricMean[node.loc,]==0]
			cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
			effectiveNumSamples[node.loc, otu.tab[, i]!=0] <- effectiveNumSamples[node.loc, otu.tab[, i]!=0] - 1
			clrcum[node.loc,cum[node.loc, ]!=0 & otu.tab[, i]!=0] <- log2(cum[node.loc, cum[node.loc, ]!=0 & otu.tab[, i]!=0]) - geometricMean[node.loc,cum[node.loc, ]!=0 & otu.tab[, i]!=0]
			node <- edge[node.loc, 1]
			node.loc <- which(edge2 == node)
		}
	}
	
	
	for (i in 2:n) {
		for (j in 1:(i-1)) {
			cum1 <- clrcum[, i]
			cum2 <- clrcum[, j]
			ind <- (cum1 + cum2) != 0
			cum1 <- cum1[ind]
			cum2 <- cum2[ind]		
			br.len2 <- br.len[ind]			
			#aitchison distance, like log fold change
			diff <- abs(cum1-cum2)	
			
			# Generalized UniFrac distance
			for(k in 1:length(alpha)){
				w <- br.len2 * (cum1 + cum2)^alpha[k]
				unifracs[i, j, k] <- unifracs[j, i, k] <- sum(diff * w) / sum(w)
			}			
						
			#	Unweighted UniFrac Distance
			cum1 <- (cum1 != 0)
			cum2 <- (cum2 != 0)			
			unifracs[i, j, (k + 1)] <- unifracs[j, i, (k + 1)] <- 
					sum(abs(cum1 - cum2) / (cum1 + cum2) * br.len2) / sum(br.len2)

		}
	}
	return(list(unifracs=unifracs))
}

 
PermanovaG <- function(formula, dat=NULL, ...) {
	# Distance based statistical test by combining multiple distance
	# matrices based on PERMANOVA procedure by taking the maximum of
	# pseudo-F statistics
	#	
	# Args:
	#		formula: left side of the formula is a three dimensional array
	#						of the supplied distance matrices as produced by GUniFrac
	#		dat: data.frame containing the covariates
	#		...: Parameter passing to adonis function
	#
	# Returns:
	# 	aov.tab (data.frame): columns - F.model, p.value
	#						 							rows - Covariates
	#
	
	save.seed <- get(".Random.seed", .GlobalEnv)
	lhs <- formula[[2]]
	lhs <- eval(lhs, dat, parent.frame())
	rhs <- as.character(formula)[3]
	f.perms <- -Inf
	f.stat <- -Inf
	for(i in 1:(dim(lhs)[3])){
		# .Random.seed is a global variable. 
		assign(".Random.seed", save.seed, .GlobalEnv)	
		Y <- as.dist(lhs[, , i])
		formula2 <- as.formula(paste("Y", "~", rhs))		
		obj <- adonis(formula2, dat, ...)		
		f.perms <- ifelse(f.perms > obj$f.perms, f.perms, obj$f.perms)		
		temp <- obj$aov.tab[1:ncol(f.perms), "F.Model"]
		f.stat <- ifelse(f.stat > temp, f.stat, temp)		
	}	
	pv <- (rowSums(t(f.perms) >= f.stat) + 1) / (nrow(f.perms) + 1)	
	aov.tab <- data.frame(F.Model = f.stat, p.value = pv )
	rownames(aov.tab) <- rownames(obj$aov.tab)[1:ncol(f.perms)]
	
	list(aov.tab=aov.tab)
}


Rarefy <- function (otu.tab, depth = min(rowSums(otu.tab))){
	# Rarefaction function: downsample to equal depth
	#	
	# Args:
	#		otu.tab: OTU count table, row - n sample, column - q OTU
	#		depth: required sequencing depth 
	#
	# Returns:
	# 	otu.tab.rff: Rarefied OTU table
	#		discard: labels of discarded samples
	#
	otu.tab <- as.matrix(otu.tab)
	ind <- (rowSums(otu.tab) < depth)
	sam.discard <- rownames(otu.tab)[ind]
	otu.tab <- otu.tab[!ind, ]
	
	rarefy <- function(x, depth){
		y <- sample(rep(1:length(x), x), depth)
		y.tab <- table(y)
		z <- numeric(length(x))
		z[as.numeric(names(y.tab))] <- y.tab
		z
	}
	otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
	rownames(otu.tab.rff) <- rownames(otu.tab)
	colnames(otu.tab.rff) <- colnames(otu.tab)
	return(list(otu.tab.rff=otu.tab.rff, discard=sam.discard))
}


