#GETTING THE DATA FROM Human Microbiome Project OTU table
# output tables for low, medium, and high sequencing depth
#samples may overlap between tables (want extra samples per group to perform replicates)

# table comes out with samples in columns and OTUs in rows

options(error=recover)

library(phangorn)
library(vegan)

writeFile <- function(saliva,stool,tree,filename,condition1,condition2,otuIDs) {
	#concatenate	
	data <- data.frame(saliva,stool)
	colnames(data) <- sub("^X", "", colnames(data))
	rownames(data) <- otuIDs

	#make condition vector
	groups <- as.factor(c(rep(condition1,ncol(saliva)),rep(condition2,ncol(stool))))

	#get rid of zero sum rows
	data.sum <- apply(data,1,sum)
	data.0 <- data[data.sum > 0,]
	data <- data.0

	# get rid of extra OTUs in tree
	tree$tip.label <- gsub("'","",tree$tip.label)
	absent <- tree$tip.label[!(tree$tip.label %in% rownames(data))]
	if (length(absent) != 0) {
			tree <- drop.tip(tree, absent)
	}
	write.tree(tree,file=paste(filename,"subtree.tre",sep="_"))

	#make sample names that contain condition

	colnames(data) <- paste(groups,colnames(data),sep="_")

	#write otu counts into table
	write.table(data,file=paste(filename,"hmp_data.txt",sep="_"),sep="\t",quote=FALSE)
	# read in with read.table("hmp_mouth_data.txt",sep="\t",header=TRUE,row.names=1)	
}



hmpData <- "../../../../fodor/otu_table_psn_v35.txt"
hmpMetadata <- "../../../../fodor/v35_map_uniquebyPSN.txt"
treeFile <- "../../../../fodor/rep_set_v35.tre"

id <- read.table(hmpMetadata, header=TRUE, sep="\t", row.names=1)
otu <- t( read.table(hmpData, header=T, sep="\t", row.names=1, check.names=FALSE) )
tree <- read.tree(treeFile)
# BODY SITES
#  [1] "Anterior_nares"               "Attached_Keratinized_gingiva"
#  [3] "Buccal_mucosa"                "Hard_palate"                 
#  [5] "Left_Antecubital_fossa"       "Left_Retroauricular_crease"  
#  [7] "Mid_vagina"                   "Palatine_Tonsils"            
#  [9] "Posterior_fornix"             "Right_Antecubital_fossa"     
# [11] "Right_Retroauricular_crease"  "Saliva"                      
# [13] "Stool"                        "Subgingival_plaque"          
# [15] "Supragingival_plaque"         "Throat"                      
# [17] "Tongue_dorsum"                "Vaginal_introitus"  

saliva <- "Saliva"
stool <- "Stool"

saliva.id <- rownames(id)[which(id$HMPbodysubsite==saliva)]
stool.id <- rownames(id)[which(id$HMPbodysubsite==stool)]

saliva.otu <- site <- otu[rownames(otu) %in% saliva.id,]
stool.otu <- site <- otu[rownames(otu) %in% stool.id,]

otuIDs <- colnames(site)


saliva.otu <- apply(saliva.otu, 1, function(x){as.numeric(x)})
stool.otu <- apply(stool.otu, 1, function(x){as.numeric(x)})

# remove all samples with read count lower than 10,000
saliva.sum <- apply(saliva.otu,2,sum)
stool.sum <- apply(stool.otu,2,sum)

saliva.low <- saliva.otu[,which(saliva.sum<3000)]
stool.low <- stool.otu[,which(stool.sum<3000)]

saliva.med <- saliva.otu[,which((saliva.sum>3000) & (saliva.sum < 6000))]
stool.med <- stool.otu[,which((stool>3000) & (stool.sum < 6000))]

saliva.high <- saliva.otu[,which(saliva.sum>6000)]
stool.high <- stool.otu[,which(stool.sum>6000)]

writeFile(saliva.low,stool.low,tree,"low_sequencing_depth",saliva,stool,otuIDs)
writeFile(saliva.med,stool.med,tree,"med_sequencing_depth",saliva,stool,otuIDs)
writeFile(saliva.high,stool.high,tree,"high_sequencing_depth",saliva,stool,otuIDs)

#examine diversity
saliva.low.div <- diversity(saliva.low)
stool.low.div <- diversity(stool.low)

saliva.med.div <- diversity(saliva.med)
stool.med.div <- diversity(stool.med)

saliva.high.div <- diversity(saliva.high)
stool.high.div <- diversity(stool.high)

summary(saliva.low.div)
summary(stool.low.div)
summary(saliva.med.div)
summary(stool.med.div)
summary(saliva.high.div)
summary(stool.high.div)
# #pick 20 random samples from each category

# bm.rand <- bm.otu[,as.integer(sample(seq(1,length(colnames(bm.otu)),1),20,replace=FALSE))]
# td.rand <- td.otu[,as.integer(sample(seq(1,length(colnames(td.otu)),1),20,replace=FALSE))]
# akg.rand <- akg.otu[,as.integer(sample(seq(1,length(colnames(akg.otu)),1),20,replace=FALSE))]
# hp.rand <- hp.otu[,as.integer(sample(seq(1,length(colnames(hp.otu)),1),20,replace=FALSE))]
# s.rand <- s.otu[,as.integer(sample(seq(1,length(colnames(s.otu)),1),20,replace=FALSE))]


