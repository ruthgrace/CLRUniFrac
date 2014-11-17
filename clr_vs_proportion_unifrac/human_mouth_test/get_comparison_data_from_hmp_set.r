#GETTING THE DATA FROM Human Microbiome Project OTU table

hmpData <- "../../../../fodor/otu_table_psn_v35.txt"
hmpMetadata <- "../../../../fodor/v35_map_uniquebyPSN.txt"
id <- read.table(hmpMetadata, header=TRUE, sep="\t", row.names=1)
otu <- t( read.table(hmpData, header=T, sep="\t", row.names=1, check.names=FALSE) )

rownames(otu) <- sub("^X", "", rownames(otu))

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

#i want to use "Buccal_mucosa" (lining of mouth), "Tongue_dorsum", "Attached_Keratinized_gingiva" (where hard and soft gum meet), "Hard_palate" , "Saliva"

#can also try "Palatine_Tonsils"  vs "Throat" later

bm <- "Buccal_mucosa"
td <- "Tongue_dorsum"
akg <- "Attached_Keratinized_gingiva"
hp <- "Hard_palate"
s <- "Saliva"

bm.id <- rownames(id)[which(id$HMPbodysubsite==bm)]
td.id <- rownames(id)[which(id$HMPbodysubsite==td)]
akg.id <- rownames(id)[which(id$HMPbodysubsite==akg)]
hp.id <- rownames(id)[which(id$HMPbodysubsite==hp)]
s.id <- rownames(id)[which(id$HMPbodysubsite==s)]

bm.otu <- site <- otu[rownames(otu) %in% bm.id,]
td.otu <- site <- otu[rownames(otu) %in% td.id,]
akg.otu <- site <- otu[rownames(otu) %in% akg.id,]
hp.otu <- site <- otu[rownames(otu) %in% hp.id,]
s.otu <- site <- otu[rownames(otu) %in% s.id,]

bm.otu <- apply(bm.otu, 1, function(x){as.numeric(x)})
td.otu <- apply(td.otu, 1, function(x){as.numeric(x)})
akg.otu <- apply(akg.otu, 1, function(x){as.numeric(x)})
hp.otu <- apply(hp.otu, 1, function(x){as.numeric(x)})
s.otu <- apply(s.otu, 1, function(x){as.numeric(x)})

#pick 20 random samples from each category

bm.rand <- bm.otu[,as.integer(sample(seq(1,length(colnames(bm.otu)),1),20,replace=FALSE))]
td.rand <- td.otu[,as.integer(sample(seq(1,length(colnames(td.otu)),1),20,replace=FALSE))]
akg.rand <- akg.otu[,as.integer(sample(seq(1,length(colnames(akg.otu)),1),20,replace=FALSE))]
hp.rand <- hp.otu[,as.integer(sample(seq(1,length(colnames(hp.otu)),1),20,replace=FALSE))]
s.rand <- s.otu[,as.integer(sample(seq(1,length(colnames(s.otu)),1),20,replace=FALSE))]

#concatenate

data <- data.frame(bm.rand,td.rand,akg.rand,hp.rand,s.rand)
colnames(data) <- sub("^X", "", colnames(data))

#make condition vector

groups <- as.factor(c(rep("bm",20),rep("td",20),rep("akg",20),rep("hp",20),rep("s",20)))

#get rid of zero sum rows

data.sum <- apply(data,1,sum)
data.0 <- data[data.sum > 0,]
data <- data.0

#make sample names that contain condition

colnames(data) <- paste(groups,colnames(data),sep="_")

#write otu counts into table

write.table(data,file="hmp_mouth_data.txt",sep="\t",quote=FALSE)

# read in with read.table("hmp_mouth_data.txt",sep="\t",header=TRUE,row.names=1)