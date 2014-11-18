
library(ape)
library(phangorn)

mouth.otu <- read.table("hmp_mouth_data.txt",sep="\t",header=TRUE,row.names=1)

groups <- as.factor(c(rep("bm",20),rep("td",20),rep("akg",20),rep("hp",20),rep("s",20)))

# read and root tree (rooted tree is required)
mouth.tree <- read.tree("./rep_set_v35.tre")
mouth.tree <- midpoint(mouth.tree)
#get rid of extra quotes on OTU labels
mouth.tree$tip.label <- gsub("'","",mouth.tree$tip.label)

source("../../CLRUniFrac.R")
source("../../GUniFrac.R")

mouth.otu <- t(mouth.otu)

clrUnifrac <- CLRUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]
gUnifrac <- GUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]

