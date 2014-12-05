options(error=recover)

library(ape)
library(phangorn)


#all CLR DIRICHLET commented out while the method is being fixed.
#attempt at sparsity filter instead -- 30 counts minimum is what is stable


mouth.otu <- read.table("hmp_mouth_data_high_read_count.txt",sep="\t",header=TRUE,row.names=1)

groups <- as.factor(c(rep("buccal mucosa",20),rep("tongue dorsum",20),rep("attached keratinized gingiva",20),rep("hard palate",20),rep("saliva",20)))

# read and root tree (rooted tree is required)
mouth.tree <- read.tree("./rep_set_v35_subtree.tre")
if (!is.rooted(mouth.tree)) {
	mouth.tree <- midpoint(mouth.tree)
}
#get rid of extra quotes on OTU labels
mouth.tree$tip.label <- gsub("'","",mouth.tree$tip.label)

#get rid of extra OTUs in tree
# absent <- mouth.tree$tip.label[!(mouth.tree$tip.label %in% colnames(mouth.otu))]
# if (length(absent) != 0) {
# 		mouth.tree <- drop.tip(mouth.tree, absent)
# 		write.tree(mouth.tree,file="rep_set_v35_subtree.tre")
# }

source("../../CLRUniFrac.R")
source("../../GUniFrac.R")
source("../../EntropyUniFrac.R")
#source("../../CLRDirichletUniFrac.R")

#format otu table for input into unifrac methods (rownames are samples, colnames are OTUs)
mouth.otu <- t(mouth.otu)

#get rid of the three OTUs that aren't in the tree (only 12 reads discarded in total)
mouth.otu <- mouth.otu[,which(colnames(mouth.otu) %in% mouth.tree$tip.label)]

#RAREFY
mouth.original <- mouth.otu
# rarefiedData <- Rarefy(mouth.otu,depth=2000)
# mouth.otu <- rarefiedData[[1]]
# mouth.original <- mouth.original[match(rownames(mouth.otu),rownames(mouth.original)),]


#order otus by abundance (least to most)
taxaOrder <- rev(order(apply(mouth.otu,2,sum)))
mouth.otu <- mouth.otu[,taxaOrder]

#sparsity filter
#remove all OTUs for which the minimum count is < 30
# mouth.otu.min <- apply(mouth.otu,2,min)
# mouth.otu <- mouth.otu[,which(mouth.otu.min) >= 30)]
# mouth.original <- mouth.original[,which(mouth.otu.min) >= 30)]

clrUnifrac <- CLRUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]
gUnifrac <- GUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]
eUnifrac <- InformationUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]
#clrDirichletUnifrac <- CLRDirichletUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]

#calculate principle coordinates of analysis
clrUnifrac.pcoa <- pcoa(clrUnifrac)
gUnifrac.pcoa <- pcoa(gUnifrac)
eUnifrac.pcoa <- pcoa(eUnifrac)
#clrDirichletUnifrac.pcoa <- pcoa(clrDirichletUnifrac)

# calculate total variance explained
clrUnifrac.varExplained <- sum(apply(clrUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
gUnifrac.varExplained <- sum(apply(gUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
eUnifrac.varExplained <- sum(apply(eUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
#clrDirichletUnifrac.varExplained <- sum(apply(clrDirichletUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))

# calculate proportion of variance explained by first component
clrUnifrac.pc1.varEx <- sd(clrUnifrac.pcoa$vector[,1])*sd(clrUnifrac.pcoa$vector[,1])/clrUnifrac.varExplained
#calculate proportion of variance explained by second component
clrUnifrac.pc2.varEx <- sd(clrUnifrac.pcoa$vector[,2])*sd(clrUnifrac.pcoa$vector[,2])/clrUnifrac.varExplained

gUnifrac.pc1.varEx <- sd(gUnifrac.pcoa$vector[,1])*sd(gUnifrac.pcoa$vector[,1])/gUnifrac.varExplained
gUnifrac.pc2.varEx <- sd(gUnifrac.pcoa$vector[,2])*sd(gUnifrac.pcoa$vector[,2])/gUnifrac.varExplained

eUnifrac.pc1.varEx <- sd(eUnifrac.pcoa$vector[,1])*sd(eUnifrac.pcoa$vector[,1])/eUnifrac.varExplained
eUnifrac.pc2.varEx <- sd(eUnifrac.pcoa$vector[,2])*sd(eUnifrac.pcoa$vector[,2])/eUnifrac.varExplained

#clrDirichletUnifrac.pc1.varEx <- sd(clrDirichletUnifrac.pcoa$vector[,1])*sd(clrDirichletUnifrac.pcoa$vector[,1])/clrDirichletUnifrac.varExplained
#clrDirichletUnifrac.pc2.varEx <- sd(clrDirichletUnifrac.pcoa$vector[,2])*sd(clrDirichletUnifrac.pcoa$vector[,2])/clrDirichletUnifrac.varExplained

#get otu proportions for barplot
mouth.prop <- t(apply(mouth.otu,1,function(x) x/sum(x)))

#get otu total read counts
mouth.sum <- apply(mouth.original,1,sum)

#convert to dist structure
clrUnifrac.dist <- as.dist(clrUnifrac)
gUnifrac.dist <- as.dist(gUnifrac)
eUnifrac.dist <- as.dist(eUnifrac)
#clrDirichletUnifrac.dist <- as.dist(clrDirichletUnifrac)

#"average" is most similar to UPGMA, apparently
clrUnifrac.dendo <- hclust(clrUnifrac.dist, method="average")
gUnifrac.dendo <- hclust(gUnifrac.dist, method="average")
eUnifrac.dendo <- hclust(eUnifrac.dist, method="average")
#clrDirichletUnifrac.dendo <- hclust(clrDirichletUnifrac.dist, method="average")


# test overlap & read count correlations
source("../metrics.r")

#overlap <- getOverlap(mouth.otu)
overlap <- vegdist(mouth.original,method="bray")
avg <- averageReadCount(mouth.original)


#put metrics matricies into single dimensional vectors for plotting
overlap.vector <- unlist(overlap[lower.tri(overlap,diag=TRUE)])
avg.vector <- unlist(avg[lower.tri(avg,diag=TRUE)])

#put distance matrices into single dimensional vectors for plotting
clrUnifrac.vector <- unlist(clrUnifrac[lower.tri(clrUnifrac,diag=TRUE)])
gUnifrac.vector <- unlist(gUnifrac[lower.tri(gUnifrac,diag=TRUE)])
eUnifrac.vector <- unlist(eUnifrac[lower.tri(eUnifrac,diag=TRUE)])
#clrDirichletUnifrac.vector <- unlist(clrDirichletUnifrac[lower.tri(clrDirichletUnifrac,diag=TRUE)])

#convert to dist structure
clrUnifrac.dist <- as.dist(clrUnifrac)
gUnifrac.dist <- as.dist(gUnifrac)
eUnifrac.dist <- as.dist(eUnifrac)
#clrDirichletUnifrac.dist <- as.dist(clrDirichletUnifrac)

#"average" is most similar to UPGMA, apparently
clrUnifrac.dendo <- hclust(clrUnifrac.dist, method="average")
gUnifrac.dendo <- hclust(gUnifrac.dist, method="average")
eUnifrac.dendo <- hclust(eUnifrac.dist, method="average")
#clrDirichletUnifrac.dendo <- hclust(clrDirichletUnifrac.dist, method="average")


# get default par
plotParameters <- par()

originalPalette <- palette()



#save to pdf
pdf("hmp_mouth_comparison_pcoa_high_read_count_no_bar_plots.pdf")

#plot overlap vs clrunifrac distance
plot(clrUnifrac.vector,overlap.vector,main="clr combination weights vs overlap")
#plot lowess line of best fit
lines(lowess(clrUnifrac.vector,overlap.vector), col="yellow") # lowess line (x,y)

#repeat for gunifrac
plot(gUnifrac.vector,overlap.vector,main="gunifrac vs overlap")
lines(lowess(gUnifrac.vector,overlap.vector), col="yellow") # lowess line (x,y)

plot(eUnifrac.vector,overlap.vector,main="eunifrac vs overlap")
lines(lowess(eUnifrac.vector,overlap.vector), col="yellow") # lowess line (x,y)

#repeat for clr dirichlet
#plot(clrDirichletUnifrac.vector,overlap.vector,main="clr dirichlet vs overlap")
#lines(lowess(clrDirichletUnifrac.vector,overlap.vector), col="yellow") # lowess line (x,y)


#plot number of reads vs unifrac distances (checking for read count bias)
plot(clrUnifrac.vector,avg.vector,main="CLR transform weighted\nUniFrac vs. sequencing depth",xlab="UniFrac distance",ylab="Average Total Read Count",col="palegreen",cex.lab=1.4,cex.main=2)
abline(fit <- lm(avg.vector ~ clrUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)
#lines(lowess(clrUnifrac.vector,avg.vector), col="darkorchid4") # lowess line (x,y)

plot(gUnifrac.vector,avg.vector,main="Proportional abundance weighted\nUniFrac vs. sequencing depth",xlab="UniFrac distance",ylab="Average Total Read Count",col="palegreen",cex.lab=1.4,cex.main=2)
abline(fit <- lm(avg.vector ~ gUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)
#lines(lowess(gUnifrac.vector,avg.vector), col="darkorchid4") # lowess line (x,y)

plot(eUnifrac.vector,avg.vector,main="Entropy weighted\nUniFrac vs. sequencing depth",xlab="UniFrac distance",ylab="Average Total Read Count",col="palegreen",cex.lab=1.4,cex.main=2)
abline(fit <- lm(avg.vector ~ eUnifrac.vector),col="darkorchid4")
print("clr vs overlap")
print(summary(fit)$r.squared)

#plot(clrDirichletUnifrac.vector,avg.vector,main="clr dirichlet vs avg")
#lines(lowess(clrDirichletUnifrac.vector,avg.vector), col="yellow") # lowess line (x,y)


plot(gUnifrac.vector,clrUnifrac.vector,main="gunifrac vs clrunifrac")
plot(gUnifrac.vector,eUnifrac.vector,main="gunifrac vs eunifrac")



#plot dendogram with bar plots
colors <- c("steelblue3","skyblue1", "indianred1", "mediumpurple1", "olivedrab3", "pink", "#FFED6F", "mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "bisque", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED")
palette(colors)

# #GG legacy code. Fix size, margins, position
# par(mfrow=c(2,1), mar=c(1, 3, 2, 1) + 0.1)
# plot(clrUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(mouth.prop[clrUnifrac.dendo$order,]), space=0,col=colors, las=2)

# plot(gUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(mouth.prop[gUnifrac.dendo$order,]), space=0,col=colors, las=2)

# plot(eUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(mouth.prop[eUnifrac.dendo$order,]), space=0,col=colors, las=2)

#plot(clrDirichletUnifrac.dendo, axes=F, ylab=NULL, ann=F)
#order the barplot 
#barplot(t(mouth.prop[clrDirichletUnifrac.dendo$order,]), space=0,col=colors, las=2)

par(plotParameters)





colors <- c("pink","red","purple","blue","orange","black")
palette(colors)
par(xpd=TRUE)
#plot pcoa plots with legend
plot(clrUnifrac.pcoa$vectors[,1],clrUnifrac.pcoa$vectors[,2], type="p",col=groups,main="CLR transform weighted UniFrac\npricipal coordinates analysis",xlab=paste("First Component", round(clrUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(clrUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(1.0,0.3,levels(groups),col=colors,pch=19)

plot(gUnifrac.pcoa$vectors[,1],gUnifrac.pcoa$vectors[,2], col=groups,main="Proportional abundance weighted UniFrac\npricipal coordinates analysis",xlab=paste("First Component", round(gUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(gUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.4,0.2,levels(groups),col=colors,pch=19)

plot(eUnifrac.pcoa$vectors[,1],eUnifrac.pcoa$vectors[,2], col=groups,main="Entropy weighted UniFrac\npricipal coordinates analysis",xlab=paste("First Component", round(eUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.4,0.2,levels(groups),col=colors,pch=19)



#plot(clrDirichletUnifrac.pcoa$vectors[,1],clrDirichletUnifrac.pcoa$vectors[,2], col=groups,main="clr dirichlet",xlab=paste("First Component", clrDirichletUnifrac.pc1.varEx,"variance explained"),ylab=paste("Second Component", clrDirichletUnifrac.pc2.varEx,"variance explained"),pch=19)
#legend(0.2,0.45,levels(groups),col=colors,pch=19)

par(plotParameters)
#plot pcoa first component vs. read count
plot(clrUnifrac.pcoa$vectors[,1],mouth.sum,main="clr combination weights vs first pcoa")
lines(lowess(clrUnifrac.pcoa$vectors[,1],mouth.sum), col="yellow") # lowess line (x,y)

plot(gUnifrac.pcoa$vectors[,1],mouth.sum,main="gunifrac vs first pcoa")
lines(lowess(gUnifrac.pcoa$vectors[,1],mouth.sum), col="yellow") # lowess line (x,y)

plot(eUnifrac.pcoa$vectors[,1],mouth.sum,main="eunifrac vs first pcoa")
lines(lowess(eUnifrac.pcoa$vectors[,1],mouth.sum), col="yellow") # lowess line (x,y)

#plot(clrDirichletUnifrac.pcoa$vectors[,1],mouth.sum,main="clr dirichlet vs avg")
#lines(lowess(clrDirichletUnifrac.pcoa$vectors[,1],mouth.sum), col="yellow") # lowess line (x,y)


#get shannon diversity average matrices
diversity <- getAvgShannonDiversity(mouth.original)
diversity.vector <- unlist(diversity[lower.tri(diversity,diag=TRUE)])

#get shannon diversity difference matrices
diversity.diff <- getShannonDiversityDiffMat(mouth.original)
#put into single dimensional vector for plotting
diversity.diff.vector <- unlist(diversity.diff[lower.tri(diversity.diff,diag=TRUE)])


darkorchid <- col2rgb("darkorchid4")
transparentdarkorchid <- rgb(darkorchid[1]/255,darkorchid[2]/255,darkorchid[3]/255,0.1)
#plot unifrac vs. shannon diversity distance matrix
plot(clrUnifrac.vector,diversity.diff.vector,main="clrunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)
plot(gUnifrac.vector,diversity.diff.vector,main="gunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.diff.vector,main="eunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)

plot(clrUnifrac.vector,diversity.vector,main="clrunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)
plot(gUnifrac.vector,diversity.vector,main="gunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.vector,main="eunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)

par(plotParameters)



dev.off()

for (i in 1:length(levels(groups))) {
	if (i<length(levels(groups))) {
		iplus1 <- i+1
		for (j in iplus1:length(levels(groups))) {
			printSeparation(clrUnifrac.pcoa,gUnifrac.pcoa,eUnifrac.pcoa,levels(groups)[i],levels(groups)[j],groups)
		}
	}	
}


palette(originalPalette)