options(error=recover)
#pcoa test script

library(ape)
library(phangorn)

# commenting out incorrect clr dirichlet

# get default par
plotParameters <- par()

#I hacked the GUniFrac script to use CLR weighting instead of proportional abundance rating
source("../../CLRUniFrac.R")

#I hacked the GUniFrac script to use CLR weighting instead of proportional abundance rating
#source("../../CLRDirichletUniFrac.R")

#this GUniFrac script was ripped straight from the GUniFrac package, with no changes.
source("../../GUniFrac.R")

source("../../EntropyUniFrac.R")


# read OTU table and format appropriately for input into UniFrac methods
brazil.otu.tab <- read.table("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- brazil.otu.tab$taxonomy
brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(brazil.otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
brazil.otu.tab <- brazil.otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
brazil.tree <- read.tree("./brazil_study_data/fasttree_all_seed_OTUs.tre")
brazil.tree <- midpoint(brazil.tree)

# read metadata
MyMeta<- read.table("./brazil_study_data/metadata_BVsamplesonly.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

# filter OTU table and metadata so that only samples which appear in both are retained
otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
otu_indicies <- otu_indicies[!is.na(otu_indicies)]
brazil.otu.tab <- brazil.otu.tab[otu_indicies,]
MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]

#run CLRUniFrac and GUniFrac for comparison, puts distance matrix in ruthClrUnifrac and gUnifrac
ruthClrUnifrac <- CLRUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
gUnifrac <- GUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
eUnifrac <- InformationUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
#clrDirichletUniFrac <- CLRDirichletUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]

#conditions (bv - bacterial vaginosis as scored by nugent/amsel, i - intermediate, n - normal/healthy)
groups <- MyMetaOrdered$n_status #levels bv, i, n
originalgroups <- groups
# change conditions so that samples which are more than 50% one taxa are colored by that taxa
otuSum <- apply(brazil.otu.tab,1,sum)
otuMax <- apply(brazil.otu.tab,1,max)
otuWhichMax <- apply(brazil.otu.tab,1,which.max)
otuDominated <- which(otuMax > otuSum/2)


otuMaxTax <- taxonomy[otuWhichMax]
otuDominated <- c(otuDominated[which(as.numeric(otuMaxTax[otuDominated])==32)],otuDominated[which(as.numeric(otuMaxTax[otuDominated])==33)])

taxonomyGroups <- as.character(groups)
taxonomyGroups[otuDominated] <- as.character(otuMaxTax[otuDominated])
taxonomyGroups <- as.factor(taxonomyGroups)

groups <- taxonomyGroups

# assign appropriate names to single taxa dominated groups
newLevels <- levels(taxonomyGroups)
splittaxa <- strsplit(levels(taxonomyGroups),split=";")

for (i in 1:length(splittaxa)) {
	if (length(splittaxa[[i]])>1) {
		newLevels[i] <- paste("L.",splittaxa[[i]][length(splittaxa[[i]])])
	}
	else {
		newLevels[i] <- splittaxa[[i]][1]
	}
}

levels(taxonomyGroups) <- newLevels


# caculate pcoa vectors
ruthClrUnifrac.pcoa <- pcoa(ruthClrUnifrac)
gUnifrac.pcoa <- pcoa(gUnifrac)
eUnifrac.pcoa <- pcoa(eUnifrac)
#clrDirichletUniFrac.pcoa <- pcoa(clrDirichletUniFrac)

# calculate total variance explained
ruthClrUnifrac.varExplained <- sum(apply(ruthClrUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
gUnifrac.varExplained <- sum(apply(gUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
eUnifrac.varExplained <- sum(apply(eUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
#clrDirichletUniFrac.varExplained <- sum(apply(clrDirichletUniFrac.pcoa$vector,2,function(x) sd(x)*sd(x)))

# calculate proportion of variance explained by first component
ruthClrUnifrac.pc1.varEx <- sd(ruthClrUnifrac.pcoa$vector[,1])*sd(ruthClrUnifrac.pcoa$vector[,1])/ruthClrUnifrac.varExplained
#calculate proportion of variance explained by second component
ruthClrUnifrac.pc2.varEx <- sd(ruthClrUnifrac.pcoa$vector[,2])*sd(ruthClrUnifrac.pcoa$vector[,2])/ruthClrUnifrac.varExplained

gUnifrac.pc1.varEx <- sd(gUnifrac.pcoa$vector[,1])*sd(gUnifrac.pcoa$vector[,1])/gUnifrac.varExplained
gUnifrac.pc2.varEx <- sd(gUnifrac.pcoa$vector[,2])*sd(gUnifrac.pcoa$vector[,2])/gUnifrac.varExplained

eUnifrac.pc1.varEx <- sd(eUnifrac.pcoa$vector[,1])*sd(eUnifrac.pcoa$vector[,1])/eUnifrac.varExplained
eUnifrac.pc2.varEx <- sd(eUnifrac.pcoa$vector[,2])*sd(eUnifrac.pcoa$vector[,2])/eUnifrac.varExplained

#clrDirichletUniFrac.pc1.varEx <- sd(clrDirichletUniFrac.pcoa$vector[,1])*sd(clrDirichletUniFrac.pcoa$vector[,1])/clrDirichletUniFrac.varExplained
#clrDirichletUniFrac.pc2.varEx <- sd(clrDirichletUniFrac.pcoa$vector[,2])*sd(clrDirichletUniFrac.pcoa$vector[,2])/clrDirichletUniFrac.varExplained


# test overlap & read count correlations
source("../metrics.r")

overlap <- getOverlap(brazil.otu.tab)
avg <- averageReadCount(brazil.otu.tab)


#put metrics matricies into single dimensional vectors for plotting
overlap.vector <- unlist(overlap[lower.tri(overlap,diag=TRUE)])
avg.vector <- unlist(avg[lower.tri(avg,diag=TRUE)])

#put distance matrices into single dimensional vectors for plotting
ruthClrUnifrac.vector <- unlist(ruthClrUnifrac[lower.tri(ruthClrUnifrac,diag=TRUE)])
gUnifrac.vector <- unlist(gUnifrac[lower.tri(gUnifrac,diag=TRUE)])
eUnifrac.vector <- unlist(eUnifrac[lower.tri(eUnifrac,diag=TRUE)])
#clrDirichletUniFrac.vector <- unlist(clrDirichletUniFrac[lower.tri(clrDirichletUniFrac,diag=TRUE)])

#get shannon diversity average matrices
diversity <- getAvgShannonDiversity(brazil.otu.tab)
diversity.vector <- unlist(diversity[lower.tri(diversity,diag=TRUE)])

#get shannon diversity difference matrices
diversity.diff <- getShannonDiversityDiffMat(brazil.otu.tab)
#put into single dimensional vector for plotting
diversity.diff.vector <- unlist(diversity.diff[lower.tri(diversity.diff,diag=TRUE)])

diversity.max <- getMaxShannonDiversity(brazil.otu.tab)
diversity.max.vector <- unlist(diversity.max[lower.tri(diversity.max,diag=TRUE)])

diversity.min <- getMinShannonDiversity(brazil.otu.tab)
diversity.min.vector <- unlist(diversity.min[lower.tri(diversity.min,diag=TRUE)])


#convert to dist structure
ruthClrUnifrac.dist <- as.dist(ruthClrUnifrac)
gUnifrac.dist <- as.dist(gUnifrac)
eUnifrac.dist <- as.dist(eUnifrac)
#clrDirichletUniFrac.dist <- as.dist(clrDirichletUniFrac)

#"average" is most similar to UPGMA, apparently
ruthClrUnifrac.dendo <- hclust(ruthClrUnifrac.dist, method="average")
gUnifrac.dendo <- hclust(gUnifrac.dist, method="average")
eUnifrac.dendo <- hclust(eUnifrac.dist, method="average")
#clrDirichletUniFrac.dendo <- hclust(clrDirichletUniFrac.dist, method="average")

#get otu proportions for barplot
brazil.prop <- t(apply(brazil.otu.tab,1,function(x) x/sum(x)))


#save plots as PDF
pdf("test_plots_with_brazil_study_data_no_bar_plots.pdf")

#plot dendogram with bar plots

#GG legacy code. Fix size, margins, position
par(mfrow=c(2,1), mar=c(1, 3, 2, 1) + 0.1,cex=0.3)

# plot(ruthClrUnifrac.dendo, axes=F, ylab=NULL, ann=F,hang=-1)

# #order the barplot 
# colors <- c("steelblue3","skyblue1", "indianred1", "mediumpurple1", "olivedrab3", "pink", "#FFED6F", "mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "bisque", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED")
# barplot(t(brazil.prop[ruthClrUnifrac.dendo$order,]), space=0,col=colors, las=2)

# plot(gUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(brazil.prop[gUnifrac.dendo$order,]), space=0,col=colors, las=2)

# plot(eUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(brazil.prop[eUnifrac.dendo$order,]), space=0,col=colors, las=2)

#plot(clrDirichletUniFrac.dendo, axes=F, ylab=NULL, ann=F)
#order the barplot 
#barplot(t(brazil.prop[clrDirichletUniFrac.dendo$order,]), space=0,col=colors, las=2)

par(plotParameters)


#for this data set, the colors represent
#  chocolate4: gardnerella vaginalis
#  darkolivegreen: prevotella bivia
#  cyan: lactobacillus crispatus
#  dodgerblue: lactobacillus iners
#  navy: lactobacillus gasseri(johnsonii)
#  magenta: streptococcus (unclassified)
#samples not dominated 50% or more by a single species:
#  red: bacterial vaginosis
#  orange: intermediate
#  blue: normal/healthy
palette(c("chocolate4","darkolivegreen","cyan","dodgerblue","navy","magenta","red","orange","blue","aquamarine"))


#plot overlap vs clrunifrac distance
plot(ruthClrUnifrac.vector,overlap.vector,main="CLR transform weighted\nUniFrac vs. overlap",col="palegreen",xlab="UniFrac distance",ylab="Overlap",cex.lab=1.4,cex.main=2)

#plot lowess line of best fit
#lines(lowess(ruthClrUnifrac.vector,overlap.vector), col="darkorchid4") # lowess line (x,y)
abline(fit <- lm(overlap.vector ~ ruthClrUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)
#repeat for gunifrac
plot(gUnifrac.vector,overlap.vector,main="Proportional abundance\nweighted UniFrac vs. overlap",col="palegreen",xlab="UniFrac distance",ylab="Overlap",cex.lab=1.4,cex.main=2)
abline(fit <- lm(overlap.vector ~ gUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)
#lines(lowess(gUnifrac.vector,overlap.vector), col="darkorchid4") # lowess line (x,y)

plot(eUnifrac.vector,overlap.vector,main="Entropy weighted\nUniFrac vs. overlap",col="palegreen",xlab="UniFrac distance",ylab="Overlap",cex.lab=1.4,cex.main=2)
abline(fit <- lm(overlap.vector ~ eUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)

#repeat for clr dirichlet unifrac
#plot(clrDirichletUniFrac.vector,overlap.vector,main="clr dirichlet vs overlap")
#lines(lowess(clrDirichletUniFrac.vector,overlap.vector), col="yellow") # lowess line (x,y)

palette(c("cyan","dodgerblue","red","orange","blue","black"))


#plot number of reads vs unifrac distances (checking for read count bias)
plot(ruthClrUnifrac.vector,avg.vector,main="CLR transform weighted\nUniFrac vs. sequencing depth",col="palegreen",xlab="UniFrac distance",ylab="Average Total Read Count",cex.lab=1.4,cex.main=2)
#lines(lowess(ruthClrUnifrac.vector,avg.vector), col="darkorchid4") # lowess line (x,y)
abline(fit <- lm(avg.vector ~ ruthClrUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)

plot(gUnifrac.vector,avg.vector,main="Proportional abundance weighted\nUniFrac vs. sequencing depth",col="palegreen",xlab="UniFrac distance",ylab="Average Total Read Count",cex.lab=1.4,cex.main=2)
#lines(lowess(gUnifrac.vector,avg.vector), col="darkorchid4") # lowess line (x,y)
abline(fit <- lm(avg.vector ~ gUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)

plot(eUnifrac.vector,avg.vector,main="Entropy weighted\nUniFrac vs. sequencing depth",col=rgb(.1,1,.1,0.1), pch=19,xlab="UniFrac distance",ylab="Average Total Read Count",cex.lab=1.4,cex.main=2)
#lines(lowess(gUnifrac.vector,avg.vector), col="darkorchid4") # lowess line (x,y)
abline(fit <- lm(avg.vector ~ eUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)

#plot(clrDirichletUniFrac.vector,avg.vector,main="clr dirichlet vs avg")
#lines(lowess(clrDirichletUniFrac.vector,avg.vector), col="yellow") # lowess line (x,y)


#plot pcoa plots with legend
plot(ruthClrUnifrac.pcoa$vectors[,1],ruthClrUnifrac.pcoa$vectors[,2], type="p",col=groups,main="CLR transformed weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(ruthClrUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(ruthClrUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(1.0,1.5,levels(taxonomyGroups),col=palette(),pch=19)

plot(gUnifrac.pcoa$vectors[,1],gUnifrac.pcoa$vectors[,2], col=groups,main="Proportional abundance weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(gUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(gUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.2,0.3,levels(taxonomyGroups),col=palette(),pch=19)

plot(eUnifrac.pcoa$vectors[,1],eUnifrac.pcoa$vectors[,2], col=groups,main="Entropy weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(eUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.2,0.3,levels(taxonomyGroups),col=palette(),pch=19)

#plot(clrDirichletUniFrac.pcoa$vectors[,1],clrDirichletUniFrac.pcoa$vectors[,2], col=groups,main="clr dirichlet unifrac",xlab=paste("First Component", clrDirichletUniFrac.pc1.varEx,"variance explained"),ylab=paste("Second Component", clrDirichletUniFrac.pc2.varEx,"variance explained"))
#legend(0.1,0.3,levels(taxonomyGroups),col=palette(),pch=1)


plot(gUnifrac.vector,ruthClrUnifrac.vector,main="gunifrac vs clrunifrac")
plot(gUnifrac.vector,eUnifrac.vector,main="gunifrac vs eunifrac")


#change color palette for qiime output data (no otu count information for dominant taxa)
# red for bacterial vaginosis, orange for intermediate, blue for normal/healthy
palette(c("red","orange","blue","black"))

# read in unifrac distances from qiime
unifracWeights <- read.table("./brazil_study_data/weighted_unifrac_dm_from_qiime.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
#match up with metadata
indicies <- match(rownames(MyMeta),rownames(unifracWeights))
indicies <- indicies[!is.na(indicies)]
unifracWeightsFiltered <- unifracWeights[indicies,indicies]
meta_indicies <- match(colnames(unifracWeightsFiltered),rownames(MyMeta))
meta_indicies <- meta_indicies[!is.na(meta_indicies)]
#set condition (bv/i/n)
groups <- MyMeta[meta_indicies,]$n_status

#calculate pcoa
unifracWeightsFiltered.ape.pcoa <- pcoa(unifracWeightsFiltered)

#calculate proportion of variance explained by 1st and 2nd components
unifracWeightsFiltered.varExplained <- sum(apply(unifracWeightsFiltered.ape.pcoa$vector,2,function(x) sd(x)*sd(x)))
unifracWeightsFiltered.pc1.varEx <- sd(unifracWeightsFiltered.ape.pcoa$vector[,1])*sd(unifracWeightsFiltered.ape.pcoa$vector[,1])/unifracWeightsFiltered.varExplained
unifracWeightsFiltered.pc2.varEx <- sd(unifracWeightsFiltered.ape.pcoa$vector[,2])*sd(unifracWeightsFiltered.ape.pcoa$vector[,2])/unifracWeightsFiltered.varExplained

#plot qiime unifrac distances pcoa with legend
plot(unifracWeightsFiltered.ape.pcoa$vectors[,1],unifracWeightsFiltered.ape.pcoa$vectors[,2], col=groups,main="pcoa from qiime unifrac distances",xlab=paste("First Component", unifracWeightsFiltered.pc1.varEx,"variance explained"),ylab=paste("Second Component", unifracWeightsFiltered.pc2.varEx,"variance explained"))
legend(-0.6,0.4,levels(groups),col=palette(),pch=1)

#match qiime pcoa vectors with metadata
qiimePCOA <- read.table("./brazil_study_data/weighted_unifrac_pc_from_qiime.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
qiimePCOA <- qiimePCOA[-nrow(qiimePCOA)+1:-nrow(qiimePCOA),]
indicies <- match(rownames(qiimePCOA),rownames(MyMeta))
indicies <- indicies[!is.na(indicies)]
groups <- MyMeta[indicies,]$n_status
qiimePCOA <- qiimePCOA[match(rownames(MyMeta[indicies,]),rownames(qiimePCOA)),]

#calculate proportion of variance explained
qiimePCOA.varExplained <- sum(apply(qiimePCOA,2,function(x) sd(x)*sd(x)))
qiimePCOA.pc1.varEx <- sd(qiimePCOA[,1])*sd(qiimePCOA[,1])/qiimePCOA.varExplained
qiimePCOA.pc2.varEx <- sd(qiimePCOA[,2])*sd(qiimePCOA[,2])/qiimePCOA.varExplained

#plot qiime pecoa vectors with legend
plot(qiimePCOA[,1],qiimePCOA[,2], col=groups,main="qiime pcoa",xlab=paste("First Component", qiimePCOA.pc1.varEx,"variance explained"),ylab=paste("Second Component", qiimePCOA.pc2.varEx,"variance explained"))
legend(0.2,0.3,levels(groups),col=palette(),pch=1)


#plot pcoa first component vs. read count
plot(ruthClrUnifrac.pcoa$vectors[,1],otuSum,main="Proportional abundance weighted UniFrac vs. PCoA first component",col="palegreen",xlab="UniFrac distance",ylab="Average total read count")
lines(lowess(ruthClrUnifrac.pcoa$vectors[,1],otuSum), col="darkorchid4") # lowess line (x,y)

plot(gUnifrac.pcoa$vectors[,1],otuSum,main="CLR transform weighted UniFrac vs. PCoA first component",col="palegreen",xlab="UniFrac distance",ylab="Average total read count")
lines(lowess(gUnifrac.pcoa$vectors[,1],otuSum), col="darkorchid4") # lowess line (x,y)

#plot(clrDirichletUniFrac.pcoa$vectors[,1],otuSum,main="clr dirichlet vs avg")
#lines(lowess(clrDirichletUniFrac.pcoa$vectors[,1],otuSum), col="yellow") # lowess line (x,y)

printSeparation(ruthClrUnifrac.pcoa,gUnifrac.pcoa,eUnifrac.pcoa,levels(originalgroups)[1],levels(originalgroups)[3],originalgroups)

palette(c("chocolate4","darkolivegreen","cyan","dodgerblue","navy","magenta","red","orange","blue","aquamarine","darkorchid4"))

darkorchid <- col2rgb("darkorchid4")
transparentdarkorchid <- rgb(darkorchid[1]/255,darkorchid[2]/255,darkorchid[3]/255,0.1)
#plot unifrac vs. shannon diversity distance matrix
plot(ruthClrUnifrac.vector,diversity.diff.vector,main="clrunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)
plot(gUnifrac.vector,diversity.diff.vector,main="gunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.diff.vector,main="eunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)

plot(ruthClrUnifrac.vector,diversity.vector,main="clrunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)
plot(gUnifrac.vector,diversity.vector,main="gunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.vector,main="eunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)

plot(ruthClrUnifrac.vector,diversity.max.vector,main="clrunifrac vs max shannon diversity",col=transparentdarkorchid, pch=19)
plot(gUnifrac.vector,diversity.max.vector,main="gunifrac vs max shannon diversity",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.max.vector,main="eunifrac vs max shannon diversity",col=transparentdarkorchid, pch=19)

plot(ruthClrUnifrac.vector,diversity.min.vector,main="clrunifrac vs min shannon diversity",col=transparentdarkorchid, pch=19)
plot(gUnifrac.vector,diversity.min.vector,main="gunifrac vs min shannon diversity",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.min.vector,main="eunifrac vs min shannon diversity",col=transparentdarkorchid, pch=19)

par(plotParameters)

dev.off()
