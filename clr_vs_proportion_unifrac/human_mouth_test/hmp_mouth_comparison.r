
library(ape)
library(phangorn)

mouth.otu <- read.table("hmp_mouth_data.txt",sep="\t",header=TRUE,row.names=1)

groups <- as.factor(c(rep("bm",20),rep("td",20),rep("akg",20),rep("hp",20),rep("s",20)))

source("../../CLRUniFrac.R")
source("../../GUniFrac.R")

mouth.otu <- t(data)

#HOW DO I GET THE TREE

clrUnifrac <- CorrectCLRUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]
gUnifrac <- GUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]
