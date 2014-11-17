data <- read.table("hmp_mouth_data.txt",sep="\t",header=TRUE,row.names=1)

groups <- as.factor(c(rep("bm",20),rep("td",20),rep("akg",20),rep("hp",20),rep("s",20)))

