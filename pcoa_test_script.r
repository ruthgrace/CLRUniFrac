#pcoa test script

library(ape)
library(phangorn)

#I hacked the GUniFrac script to use CLR weighting instead of proportional abundance rating
source("../CLRUniFrac.R")

#this GUniFrac script was ripped straight from the GUniFrac package, with no changes.
source("../GUniFrac.R")

brazil.otu.tab <- read.table("../brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
taxonomy <- brazil.otu.tab$taxonomy
brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

brazil.tree <- read.tree("../brazil_study_data/fasttree_all_seed_OTUs.tre")
brazil.tree <- midpoint(brazil.tree)

MyMeta<- read.table("../brazil_study_data/metadata_BVsamplesonly.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
otu_indicies <- otu_indicies[!is.na(otu_indicies)]

brazil.otu.tab <- brazil.otu.tab[otu_indicies,]

MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]


ruthClrUnifrac <- CorrectCLRUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
gUnifrac <- GUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]

groups <- MyMetaOrdered$n_status #levels bv, i, n

otuSum <- apply(brazil.otu.tab,1,sum)
otuMax <- apply(brazil.otu.tab,1,max)
otuWhichMax <- apply(brazil.otu.tab,1,which.max)
otuDominated <- which(otuMax > otuSum/2)
otuMaxTax <- taxonomy[otuWhichMax]

taxonomyGroups <- as.character(groups)
taxonomyGroups[otuDominated] <- as.character(otuMaxTax[otuDominated])
taxonomyGroups <- as.factor(taxonomyGroups)

groups <- taxonomyGroups


ruthClrUnifrac.pca <- pcoa(ruthClrUnifrac)
gUnifrac.pca <- pcoa(gUnifrac)

ruthClrUnifrac.varExplained <- sum(apply(ruthClrUnifrac.pca$vector,2,function(x) sd(x)*sd(x)))
gUnifrac.varExplained <- sum(apply(gUnifrac.pca$vector,2,function(x) sd(x)*sd(x)))

ruthClrUnifrac.pc1.varEx <- sd(ruthClrUnifrac.pca$vector[,1])*sd(ruthClrUnifrac.pca$vector[,1])/ruthClrUnifrac.varExplained
ruthClrUnifrac.pc2.varEx <- sd(ruthClrUnifrac.pca$vector[,2])*sd(ruthClrUnifrac.pca$vector[,2])/ruthClrUnifrac.varExplained

gUnifrac.pc1.varEx <- sd(gUnifrac.pca$vector[,1])*sd(gUnifrac.pca$vector[,1])/gUnifrac.varExplained
gUnifrac.pc2.varEx <- sd(gUnifrac.pca$vector[,2])*sd(gUnifrac.pca$vector[,2])/gUnifrac.varExplained


# test overlap & read count correlations

source("metrics.r")

# overlap <- getOverlap(MyOTU)
# avg <- averageReadCount(MyOTU)

overlap <- getOverlap(brazil.otu.tab)

avg <- averageReadCount(brazil.otu.tab)

# overlap <- getOverlap(throat.otu.tab)

# avg <- averageReadCount(throat.otu.tab)

newLevels <- levels(taxonomyGroups)
splittaxa <- strsplit(levels(taxonomyGroups),split=";")

for (i in 1:length(splittaxa)) {
	if (length(splittaxa[[i]])>1) {
		newLevels[i] <- paste(splittaxa[[i]][length(splittaxa[[i]])-1],splittaxa[[i]][length(splittaxa[[i]])])
	}
	else {
		newLevels[i] <- splittaxa[[i]][1]
	}
}

levels(taxonomyGroups) <- newLevels

#triangle your data and put into list
overlap.vector <- unlist(overlap[lower.tri(overlap,diag=TRUE)])
avg.vector <- unlist(avg[lower.tri(avg,diag=TRUE)])

ruthClrUnifrac.vector <- unlist(ruthClrUnifrac[lower.tri(ruthClrUnifrac,diag=TRUE)])
gUnifrac.vector <- unlist(gUnifrac[lower.tri(gUnifrac,diag=TRUE)])

pdf("test_plots_with_brazil_study_data.pdf")

palette(c("chocolate4","darkolivegreen","cyan","dodgerblue","navy","magenta","red","orange","blue","black"))

#plot
plot(ruthClrUnifrac.vector,overlap.vector,main="clr combination weights vs overlap")

lines(lowess(ruthClrUnifrac.vector,overlap.vector), col="yellow") # lowess line (x,y)

plot(gUnifrac.vector,overlap.vector,main="gunifrac vs overlap")

lines(lowess(gUnifrac.vector,overlap.vector), col="yellow") # lowess line (x,y)

plot(ruthClrUnifrac.vector,avg.vector,main="clr combination weights vs avg")
plot(gUnifrac.vector,avg.vector,main="gunifrac vs avg")

plot(ruthClrUnifrac.pca$vectors[,1],ruthClrUnifrac.pca$vectors[,2], type="p",col=groups,main="clr combination weights",xlab=paste("First Component", ruthClrUnifrac.pc1.varEx,"variance explained"),ylab=paste("Second Component", ruthClrUnifrac.pc2.varEx,"variance explained"))

legend(0.5,1.5,levels(taxonomyGroups),col=palette(),pch=1)



plot(gUnifrac.pca$vectors[,1],gUnifrac.pca$vectors[,2], col=groups,main="gunifrac",xlab=paste("First Component", gUnifrac.pc1.varEx,"variance explained"),ylab=paste("Second Component", gUnifrac.pc2.varEx,"variance explained"))

legend(0.1,0.3,levels(taxonomyGroups),col=palette(),pch=1)


palette(c("red","orange","blue","black"))

unifracWeights <- read.table("./brazil_study_data/weighted_unifrac_dm_from_qiime.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
indicies <- match(rownames(MyMeta),rownames(unifracWeights))
indicies <- indicies[!is.na(indicies)]

unifracWeightsFiltered <- unifracWeights[indicies,indicies]

meta_indicies <- match(colnames(unifracWeightsFiltered),rownames(MyMeta))
meta_indicies <- meta_indicies[!is.na(meta_indicies)]
groups <- MyMeta[meta_indicies,]$n_status

unifracWeightsFiltered.ape.pcoa <- pcoa(unifracWeightsFiltered)


unifracWeightsFiltered.varExplained <- sum(apply(unifracWeightsFiltered.ape.pcoa$vector,2,function(x) sd(x)*sd(x)))

unifracWeightsFiltered.pc1.varEx <- sd(unifracWeightsFiltered.ape.pcoa$vector[,1])*sd(unifracWeightsFiltered.ape.pcoa$vector[,1])/unifracWeightsFiltered.varExplained
unifracWeightsFiltered.pc2.varEx <- sd(unifracWeightsFiltered.ape.pcoa$vector[,2])*sd(unifracWeightsFiltered.ape.pcoa$vector[,2])/unifracWeightsFiltered.varExplained

plot(unifracWeightsFiltered.ape.pcoa$vectors[,1],unifracWeightsFiltered.ape.pcoa$vectors[,2], col=groups,main="pcoa from qiime unifrac distances",xlab=paste("First Component", unifracWeightsFiltered.pc1.varEx,"variance explained"),ylab=paste("Second Component", unifracWeightsFiltered.pc2.varEx,"variance explained"))

legend(-0.6,0.4,levels(groups),col=palette(),pch=1)


qiimePCOA <- read.table("./brazil_study_data/weighted_unifrac_pc_from_qiime.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
qiimePCOA <- qiimePCOA[-nrow(qiimePCOA)+1:-nrow(qiimePCOA),]
indicies <- match(rownames(qiimePCOA),rownames(MyMeta))
indicies <- indicies[!is.na(indicies)]
groups <- MyMeta[indicies,]$n_status
qiimePCOA <- qiimePCOA[match(rownames(MyMeta[indicies,]),rownames(qiimePCOA)),]


qiimePCOA.varExplained <- sum(apply(qiimePCOA,2,function(x) sd(x)*sd(x)))

qiimePCOA.pc1.varEx <- sd(qiimePCOA[,1])*sd(qiimePCOA[,1])/qiimePCOA.varExplained
qiimePCOA.pc2.varEx <- sd(qiimePCOA[,2])*sd(qiimePCOA[,2])/qiimePCOA.varExplained


plot(qiimePCOA[,1],qiimePCOA[,2], col=groups,main="qiime pcoa",xlab=paste("First Component", qiimePCOA.pc1.varEx,"variance explained"),ylab=paste("Second Component", qiimePCOA.pc2.varEx,"variance explained"))

legend(0.2,0.3,levels(groups),col=palette(),pch=1)


dev.off()
