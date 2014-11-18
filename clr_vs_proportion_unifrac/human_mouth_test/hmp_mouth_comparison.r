
library(ape)
library(phangorn)

mouth.otu <- read.table("hmp_mouth_data.txt",sep="\t",header=TRUE,row.names=1)

groups <- as.factor(c(rep("bm",20),rep("td",20),rep("akg",20),rep("hp",20),rep("s",20)))

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

#format otu table for input into unifrac methods (rownames are samples, colnames are OTUs)
mouth.otu <- t(mouth.otu)

#get rid of the three OTUs that aren't in the tree (only 12 reads discarded in total)
mouth.otu <- mouth.otu[,which(colnames(mouth.otu) %in% mouth.tree$tip.label)]

#order otus by abundance (least to most)
taxaOrder <- rev(order(apply(mouth.otu,2,sum)))
mouth.otu <- mouth.otu[,taxaOrder]


clrUnifrac <- CLRUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]
gUnifrac <- GUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]

#calculate principle coordinates of analysis
clrUnifrac.pcoa <- pcoa(clrUnifrac)
gUnifrac.pcoa <- pcoa(gUnifrac)

# calculate total variance explained
clrUnifrac.varExplained <- sum(apply(clrUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
gUnifrac.varExplained <- sum(apply(gUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))

# calculate proportion of variance explained by first component
clrUnifrac.pc1.varEx <- sd(clrUnifrac.pcoa$vector[,1])*sd(clrUnifrac.pcoa$vector[,1])/clrUnifrac.varExplained
clrUnifrac.pc2.varEx <- sd(clrUnifrac.pcoa$vector[,2])*sd(clrUnifrac.pcoa$vector[,2])/clrUnifrac.varExplained

#calculate proportion of variance explained by second component
gUnifrac.pc1.varEx <- sd(gUnifrac.pcoa$vector[,1])*sd(gUnifrac.pcoa$vector[,1])/gUnifrac.varExplained
gUnifrac.pc2.varEx <- sd(gUnifrac.pcoa$vector[,2])*sd(gUnifrac.pcoa$vector[,2])/gUnifrac.varExplained

#get otu proportions for barplot
mouth.prop <- t(apply(mouth.otu,1,function(x) x/sum(x)))

#convert to dist structure
clrUnifrac.dist <- as.dist(clrUnifrac)
gUnifrac.dist <- as.dist(gUnifrac)

#"average" is most similar to UPGMA, apparently
clrUnifrac.dendo <- hclust(clrUnifrac.dist, method="average")
gUnifrac.dendo <- hclust(gUnifrac.dist, method="average")



# test overlap & read count correlations
source("../metrics.r")

overlap <- getOverlap(mouth.otu)
avg <- averageReadCount(mouth.otu)


#put metrics matricies into single dimensional vectors for plotting
overlap.vector <- unlist(overlap[lower.tri(overlap,diag=TRUE)])
avg.vector <- unlist(avg[lower.tri(avg,diag=TRUE)])

#put distance matrices into single dimensional vectors for plotting
clrUnifrac.vector <- unlist(clrUnifrac[lower.tri(clrUnifrac,diag=TRUE)])
gUnifrac.vector <- unlist(gUnifrac[lower.tri(gUnifrac,diag=TRUE)])

#convert to dist structure
clrUnifrac.dist <- as.dist(clrUnifrac)
gUnifrac.dist <- as.dist(gUnifrac)

#"average" is most similar to UPGMA, apparently
clrUnifrac.dendo <- hclust(clrUnifrac.dist, method="average")
gUnifrac.dendo <- hclust(gUnifrac.dist, method="average")




# get default par
plotParameters <- par()

originalPalette <- palette()



#save to pdf
pdf("hmp_mouth_comparison_pcoa_no_overlap.pdf")

#plot overlap vs clrunifrac distance
plot(clrUnifrac.vector,overlap.vector,main="clr combination weights vs overlap")
#plot lowess line of best fit
lines(lowess(clrUnifrac.vector,overlap.vector), col="yellow") # lowess line (x,y)

#repeat for gunifrac
plot(gUnifrac.vector,overlap.vector,main="gunifrac vs overlap")
lines(lowess(gUnifrac.vector,overlap.vector), col="yellow") # lowess line (x,y)

#plot number of reads vs unifrac distances (checking for read count bias)
plot(clrUnifrac.vector,avg.vector,main="clr combination weights vs avg")
plot(gUnifrac.vector,avg.vector,main="gunifrac vs avg")


#plot dendogram with bar plots
colors <- c("steelblue3","skyblue1", "indianred1", "mediumpurple1", "olivedrab3", "pink", "#FFED6F", "mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "bisque", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED")
palette(colors)

#GG legacy code. Fix size, margins, position
par(mfrow=c(2,1), mar=c(1, 3, 2, 1) + 0.1)
plot(clrUnifrac.dendo, axes=F, ylab=NULL, ann=F)
#order the barplot 
barplot(t(mouth.prop[clrUnifrac.dendo$order,]), space=0,col=colors, las=2)

par(mfrow=c(2,1), mar=c(1, 3, 2, 1) + 0.1)
plot(gUnifrac.dendo, axes=F, ylab=NULL, ann=F)
#order the barplot 
barplot(t(mouth.prop[gUnifrac.dendo$order,]), space=0,col=colors, las=2)

par(plotParameters)





colors <- c("pink","red","purple","blue","orange","black")
palette(colors)

#plot pcoa plots with legend
plot(clrUnifrac.pcoa$vectors[,1],clrUnifrac.pcoa$vectors[,2], type="p",col=groups,main="clr combination weights",xlab=paste("First Component", clrUnifrac.pc1.varEx,"variance explained"),ylab=paste("Second Component", clrUnifrac.pc2.varEx,"variance explained"),pch=19)
legend(0.6,0.5,levels(groups),col=colors,pch=19)

plot(gUnifrac.pcoa$vectors[,1],gUnifrac.pcoa$vectors[,2], col=groups,main="gunifrac",xlab=paste("First Component", gUnifrac.pc1.varEx,"variance explained"),ylab=paste("Second Component", gUnifrac.pc2.varEx,"variance explained"),pch=19)
legend(0.2,0.45,levels(groups),col=colors,pch=19)

dev.off()

palette(originalPalette)