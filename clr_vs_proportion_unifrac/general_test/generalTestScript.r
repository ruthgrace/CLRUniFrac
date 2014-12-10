options(error=recover)

otuFile <- "td_OTU_tag_mapped_lineage.txt"
treeFile <- "fasttree_all_seed_OTUs.tre"
metaDataFile <- "metadata.txt"


library(ape)
library(phangorn)


# get default par
plotParameters <- par()

source("../../GUniFrac.R")
source("../../EntropyUniFrac.R")

# read OTU table and format appropriately for input into UniFrac methods
otu.tab <- read.table(otuFile, header=T, sep="\t", row.names=1, check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- otu.tab$taxonomy
otu.tab <- otu.tab[-length(colnames(otu.tab))]
otu.tab <- t(as.matrix(otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
otu.tab <- otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
tree <- read.tree(treeFile)
tree <- midpoint(tree)

# read metadata
MyMeta<- read.table(metaDataFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

# filter OTU table and metadata so that only samples which appear in both are retained
otu_indicies <- match(rownames(MyMeta),rownames(otu.tab))
otu_indicies <- otu_indicies[!is.na(otu_indicies)]
otu.tab <- otu.tab[otu_indicies,]
MyMetaOrdered <- MyMeta[match(rownames(otu.tab),rownames(MyMeta)),]

gUnifrac <- GUniFrac(otu.tab, tree, alpha = c(1))

wUnifrac <- gUnifrac$unifrac[,,1]
uwUnifrac <- gUnifrac$unifrac[,,3]
eUnifrac <- InformationUniFrac(otu.tab, tree, alpha = c(1))$unifrac[,,1]

groups <- MyMetaOrdered$Age #levels bv, i, n
originalgroups <- groups

groups[which(groups<=17)] <- 1
groups[which(groups>17)] <- 2
groups <- as.factor(groups)

otuSum <- apply(otu.tab,1,sum)

# caculate pcoa vectors
wUnifrac.pcoa <- pcoa(wUnifrac)
uwUnifrac.pcoa <- pcoa(uwUnifrac)
eUnifrac.pcoa <- pcoa(eUnifrac)
#clrDirichletUniFrac.pcoa <- pcoa(clrDirichletUniFrac)

# calculate total variance explained
wUnifrac.varExplained <- sum(apply(wUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
uwUnifrac.varExplained <- sum(apply(uwUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
eUnifrac.varExplained <- sum(apply(eUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
#clrDirichletUniFrac.varExplained <- sum(apply(clrDirichletUniFrac.pcoa$vector,2,function(x) sd(x)*sd(x)))

# calculate proportion of variance explained by first component
wUnifrac.pc1.varEx <- sd(wUnifrac.pcoa$vector[,1])*sd(wUnifrac.pcoa$vector[,1])/wUnifrac.varExplained
#calculate proportion of variance explained by second component
wUnifrac.pc2.varEx <- sd(wUnifrac.pcoa$vector[,2])*sd(wUnifrac.pcoa$vector[,2])/wUnifrac.varExplained

uwUnifrac.pc1.varEx <- sd(uwUnifrac.pcoa$vector[,1])*sd(uwUnifrac.pcoa$vector[,1])/uwUnifrac.varExplained
uwUnifrac.pc2.varEx <- sd(uwUnifrac.pcoa$vector[,2])*sd(uwUnifrac.pcoa$vector[,2])/uwUnifrac.varExplained

eUnifrac.pc1.varEx <- sd(eUnifrac.pcoa$vector[,1])*sd(eUnifrac.pcoa$vector[,1])/eUnifrac.varExplained
eUnifrac.pc2.varEx <- sd(eUnifrac.pcoa$vector[,2])*sd(eUnifrac.pcoa$vector[,2])/eUnifrac.varExplained

# test overlap & read count correlations
source("../metrics.r")

overlap <- getOverlap(otu.tab)
avg <- averageReadCount(otu.tab)


#put metrics matricies into single dimensional vectors for plotting
overlap.vector <- unlist(overlap[lower.tri(overlap,diag=TRUE)])
avg.vector <- unlist(avg[lower.tri(avg,diag=TRUE)])

#put distance matrices into single dimensional vectors for plotting
wUnifrac.vector <- unlist(wUnifrac[lower.tri(wUnifrac,diag=TRUE)])
uwUnifrac.vector <- unlist(uwUnifrac[lower.tri(uwUnifrac,diag=TRUE)])
eUnifrac.vector <- unlist(eUnifrac[lower.tri(eUnifrac,diag=TRUE)])

#get shannon diversity average matrices
diversity <- getAvgShannonDiversity(otu.tab)
diversity.vector <- unlist(diversity[lower.tri(diversity,diag=TRUE)])

#get shannon diversity difference matrices
diversity.diff <- getShannonDiversityDiffMat(otu.tab)
#put into single dimensional vector for plotting
diversity.diff.vector <- unlist(diversity.diff[lower.tri(diversity.diff,diag=TRUE)])

diversity.max <- getMaxShannonDiversity(otu.tab)
diversity.max.vector <- unlist(diversity.max[lower.tri(diversity.max,diag=TRUE)])

diversity.min <- getMinShannonDiversity(otu.tab)
diversity.min.vector <- unlist(diversity.min[lower.tri(diversity.min,diag=TRUE)])


#convert to dist structure
wUnifrac.dist <- as.dist(wUnifrac)
uwUnifrac.dist <- as.dist(uwUnifrac)
eUnifrac.dist <- as.dist(eUnifrac)

#"average" is most similar to UPGMA, apparently
wUnifrac.dendo <- hclust(wUnifrac.dist, method="average")
uwUnifrac.dendo <- hclust(uwUnifrac.dist, method="average")
eUnifrac.dendo <- hclust(eUnifrac.dist, method="average")

#get otu proportions for barplot
otu.prop <- t(apply(otu.tab,1,function(x) x/sum(x)))

#save plots as PDF
pdf("general_test_plots.pdf")

#plot dendogram with bar plots

#GG legacy code. Fix size, margins, position
par(mfrow=c(2,1), mar=c(1, 3, 2, 1) + 0.1,cex=0.3)
colors <- c("steelblue3","skyblue1", "indianred1", "mediumpurple1", "olivedrab3", "pink", "#FFED6F", "mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "bisque", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED")

plot(wUnifrac.dendo, axes=F, ylab=NULL, ann=F,hang=-1)
#order the barplot 
barplot(t(otu.prop[wUnifrac.dendo$order,]), space=0,col=colors, las=2)

plot(uwUnifrac.dendo, axes=F, ylab=NULL, ann=F)
#order the barplot 
barplot(t(otu.prop[uwUnifrac.dendo$order,]), space=0,col=colors, las=2)

plot(eUnifrac.dendo, axes=F, ylab=NULL, ann=F)
#order the barplot 
barplot(t(otu.prop[eUnifrac.dendo$order,]), space=0,col=colors, las=2)


par(plotParameters)


#for vaginal data sets, the colors represent
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

#comparison with overlap
plot(wUnifrac.vector,overlap.vector,main="weighted UniFrac\nvs. overlap",col="palegreen",xlab="UniFrac distance",ylab="Overlap",cex.lab=1.4,cex.main=2)
abline(fit <- lm(overlap.vector ~ wUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)

plot(uwUnifrac.vector,overlap.vector,main="unweighted UniFrac\nvs. overlap",col="palegreen",xlab="UniFrac distance",ylab="Overlap",cex.lab=1.4,cex.main=2)
abline(fit <- lm(overlap.vector ~ uwUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)

plot(eUnifrac.vector,overlap.vector,main="Entropy weighted\nUniFrac vs. overlap",col="palegreen",xlab="UniFrac distance",ylab="Overlap",cex.lab=1.4,cex.main=2)
abline(fit <- lm(overlap.vector ~ eUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)

#repeat for clr dirichlet unifrac
#plot(clrDirichletUniFrac.vector,overlap.vector,main="clr dirichlet vs overlap")
#lines(lowess(clrDirichletUniFrac.vector,overlap.vector), col="yellow") # lowess line (x,y)

palette(c("cyan","dodgerblue","red","orange","blue","black"))


#plot number of reads vs unifrac distances (checking for read count bias)
plot(wUnifrac.vector,avg.vector,main="weighted UniFrac\nvs. sequencing depth",col="palegreen",xlab="UniFrac distance",ylab="Average Total Read Count",cex.lab=1.4,cex.main=2)
#lines(lowess(ruthClrUnifrac.vector,avg.vector), col="darkorchid4") # lowess line (x,y)
abline(fit <- lm(avg.vector ~ wUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)

plot(uwUnifrac.vector,avg.vector,main="unweighted UniFrac\nvs. sequencing depth",col="palegreen",xlab="UniFrac distance",ylab="Average Total Read Count",cex.lab=1.4,cex.main=2)
#lines(lowess(gUnifrac.vector,avg.vector), col="darkorchid4") # lowess line (x,y)
abline(fit <- lm(avg.vector ~ uwUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)

plot(eUnifrac.vector,avg.vector,main="Entropy weighted\nUniFrac vs. sequencing depth",col=rgb(.1,1,.1,0.1), pch=19,xlab="UniFrac distance",ylab="Average Total Read Count",cex.lab=1.4,cex.main=2)
#lines(lowess(gUnifrac.vector,avg.vector), col="darkorchid4") # lowess line (x,y)
abline(fit <- lm(avg.vector ~ eUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)



#plot pcoa plots with legend
plot(wUnifrac.pcoa$vectors[,1],wUnifrac.pcoa$vectors[,2], type="p",col=groups,main="weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(wUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(wUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(1.0,1.5,levels(groups),col=palette(),pch=19)

plot(uwUnifrac.pcoa$vectors[,1],uwUnifrac.pcoa$vectors[,2], col=groups,main="unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(uwUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(uwUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.2,0.3,levels(groups),col=palette(),pch=19)

plot(eUnifrac.pcoa$vectors[,1],eUnifrac.pcoa$vectors[,2], col=groups,main="Entropy weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(eUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.2,0.3,levels(groups),col=palette(),pch=19)



plot(wUnifrac.vector,uwUnifrac.vector,main="weighted unifrac vs unweighted unifrac")
plot(wUnifrac.vector,eUnifrac.vector,main="weighted unifrac vs eunifrac")


#plot pcoa first component vs. read count
plot(wUnifrac.pcoa$vectors[,1],otuSum,main="weighted UniFrac vs. PCoA first component",col="palegreen",xlab="UniFrac distance",ylab="Average total read count")
lines(lowess(wUnifrac.pcoa$vectors[,1],otuSum), col="darkorchid4") # lowess line (x,y)

plot(uwUnifrac.pcoa$vectors[,1],otuSum,main="unweighted UniFrac vs. PCoA first component",col="palegreen",xlab="UniFrac distance",ylab="Average total read count")
lines(lowess(uwUnifrac.pcoa$vectors[,1],otuSum), col="darkorchid4") # lowess line (x,y)

plot(eUnifrac.pcoa$vectors[,1],otuSum,main="entropy weighted UniFrac vs. PCoA first component",col="palegreen",xlab="UniFrac distance",ylab="Average total read count")
lines(lowess(eUnifrac.pcoa$vectors[,1],otuSum), col="darkorchid4") # lowess line (x,y)

palette(c("chocolate4","darkolivegreen","cyan","dodgerblue","navy","magenta","red","orange","blue","aquamarine","darkorchid4"))

darkorchid <- col2rgb("darkorchid4")
transparentdarkorchid <- rgb(darkorchid[1]/255,darkorchid[2]/255,darkorchid[3]/255,0.1)

#plot unifrac vs. shannon diversity distance matrix
plot(wUnifrac.vector,diversity.diff.vector,main="clrunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)
plot(uwUnifrac.vector,diversity.diff.vector,main="gunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.diff.vector,main="eunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)

plot(wUnifrac.vector,diversity.vector,main="clrunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)
plot(uwUnifrac.vector,diversity.vector,main="gunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.vector,main="eunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)

plot(wUnifrac.vector,diversity.max.vector,main="clrunifrac vs max shannon diversity",col=transparentdarkorchid, pch=19)
plot(uwUnifrac.vector,diversity.max.vector,main="gunifrac vs max shannon diversity",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.max.vector,main="eunifrac vs max shannon diversity",col=transparentdarkorchid, pch=19)

plot(wUnifrac.vector,diversity.min.vector,main="clrunifrac vs min shannon diversity",col=transparentdarkorchid, pch=19)
plot(uwUnifrac.vector,diversity.min.vector,main="gunifrac vs min shannon diversity",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.min.vector,main="eunifrac vs min shannon diversity",col=transparentdarkorchid, pch=19)

par(plotParameters)

dev.off()
