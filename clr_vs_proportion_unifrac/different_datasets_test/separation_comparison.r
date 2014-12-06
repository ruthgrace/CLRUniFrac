# consider testing with different rarification methods later ?

options(error=recover)

library(ape)
library(phangorn)


removeTreeTipLabelSingleQuotes <- function(tree) {
	tree$tip.label <- gsub("'","",tree$tip.label)
	return(tree)
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

#source("../../CLRUniFrac.R")
source("../../GUniFrac.R")
source("../../EntropyUniFrac.R")
#source("../../CLRDirichletUniFrac.R")

#format otu table for input into unifrac methods (rownames are samples, colnames are OTUs)
low.data.t <- t(low.data)
med.data.t <- t(med.data)
high.data.t <- t(high.data)

#get rid of any OTUs that aren't in the tree (only ___ reads discarded in total)
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

# SEQUENCING DEPTH TEST

##















#sparsity filter
#remove all OTUs for which the minimum count is < 30
# mouth.otu.min <- apply(mouth.otu,2,min)
# mouth.otu <- mouth.otu[,which(mouth.otu.min) >= 30)]
# mouth.original <- mouth.original[,which(mouth.otu.min) >= 30)]


