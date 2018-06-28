#4/3/18
#Relative abundance figure of OTUs identified by LEfSe as differentially abundant in colon, aquamin Met vs. Ctrl
#Based on Anna Seekatz's code for Rush CX manuscript figure
#6/18/18 
#Obtained from Christine Bassis, Ph.D. 
library(plyr) #tool for splitting, applying, and combining data
library(stringr) #loaded to use str_replace_all() which removes all special characters in a string
library(tm) #loaded to use removeNumbers() which removes any number in a string

setwd("~/Downloads/Work/MothurFiles/C7C9_Combined/Abundance_LEFse")

#take LDA found otus and convert shared file to relative abundances
shared <- read.table(file='combined.final.p19.shared', header=TRUE, row.names = 2) 
otu <- subset(shared, select = -c(label, numOtus))
otu.filtered <- subset(otu, select = c(Otu0008, Otu0003, Otu0020, Otu0041, Otu0017, Otu0073, Otu0033, Otu0018, Otu0037, Otu0089, Otu0032, Otu0025, Otu0048, Otu0091, Otu0058, Otu0031, Otu0072)) #otus selected from LDA values and tests

##TAXONOMY INFO.## 
##following code changes label "OTU0001" to corresponding genus/phyla.##
#get taxonomy file for data on otus
tax <- read.table("combined.final.0.03.cons.taxonomy",
                  row.names = 1,
                  header=TRUE,
                  check.names=FALSE,
                  comment.char="")

#pull just taxonomy string 
taxonomy <- as.vector(tax[,"Taxonomy"])

#define length of vector for the "for loop"
length <- as.numeric(length(taxonomy))

#changes lengthy taxonomic info to just a string of the phlya and genus
#sapply(): applys a function across a string
#strsplit(): splits string based on the character ";" and selects the 3rd value or 4th value (phlya and genus respectively)
for(i in 1:length){
  string <- as.character(taxonomy[i]) #maintains character status
  phyla <- sapply(strsplit(string,";"), `[`, 5) #should be 5
  genus <- sapply(strsplit(string,";"), `[`, 6) #should be 6
  #pre.pre.name <- str_c(phyla, " ", genus)#combines values  ###what it should be 
  pre.pre.pre.name <- str_c(genus) ###simple version for now
  pre.pre.name <- removeNumbers(pre.pre.pre.name) #remove numbers
  pre.name <- str_replace_all(pre.pre.name, "_unclassified", "")
  name <- str_replace_all(pre.name, "[[:punct:]]", "") #removes "(" ")" and all special characters
  taxonomy[i] <- name
}

#combines old taxonomy file and new taxonomy vector with right labels
taxonomy.complete <- cbind(tax, taxonomy)

#define length of full data set for the "for loop"
length.data <- as.numeric(dim(otu.filtered)[2])

#transposes taxonomy.complete so we can call the otu column titles
taxonomy.complete <- as.data.frame(t(taxonomy.complete))

#makes change to full data set.
#i is big data set, j is taxonomy file
#takes one otu in the big file and searches the taxonomy file until there is a match.  once it gets a match, it pulls the string we created in the last for loop and makes that the new column name
for(i in 1:length.data){
  for(j in 1:length){
    if(colnames(otu.filtered)[i] == colnames(taxonomy.complete)[j]){
      colnames(otu.filtered)[i] <- as.character(taxonomy.complete["taxonomy", j])
    }
    else{
    }
  }
}
##END TAXONOMY CODE##

#finish converting shared file to relative abundances
otu.matrix <- as.matrix(otu.filtered)
otu.rel <- otu.matrix/rowSums(otu.matrix)
write.table(otu.rel, "combined.final.p19.rel.txt", quote=FALSE , sep="\t", col.names=NA)

#combine design files and relative abundance files
all.otu.rel <- read.table(file="combined.final.p19.rel.txt", sep="\t", header=TRUE, row.names=1)
otu.design <- read.table(file = "combined_P19.design.txt", header = TRUE)
otu.design.otu.rel <- merge(otu.design, all.otu.rel, by.x=c("Sample"), by.y=c("row.names"))
write.table(otu.design.otu.rel, "otu.p19.design.otu.rel.txt", quote=FALSE, sep="\t", col.names=NA)

#deleted columns A and C in Excel, saved as otu.p19.design.otu.rel.txt
otu.p19.otu.rel <- read.table(file="otu.p19.design.otu.rel.txt", header=TRUE, sep="\t", row.names=2)
otu.p19.otu.rel <- subset(otu.p19.otu.rel, select = -c(X))
Met <- as.character(otu.design[otu.design$Group %in% c("M"), c("Sample")])
Ctrl <- as.character(otu.design[otu.design$Group %in% c("C"), c("Sample")])
Met.rel <- otu.p19.otu.rel[row.names(otu.p19.otu.rel) %in% Met, ]
Ctrl.rel <- otu.p19.otu.rel[row.names(otu.p19.otu.rel) %in% Ctrl, ]
Met.rel <- subset(Met.rel, select = -c(Group))
Ctrl.rel <- subset(Ctrl.rel, select = -c(Group))
# note: if you want logs for easier viewing:
log_Met.rel <- log10(Met.rel + 1)
log_Ctrl.rel <- log10(Ctrl.rel + 1)
log_Met.rel <- log_Met.rel[,order(colnames(log_Met.rel),decreasing=TRUE)]
log_Ctrl.rel <- log_Ctrl.rel[,order(colnames(log_Ctrl.rel),decreasing=TRUE)]


# step 3: plot
# OTU abundance differences: LOG transformed
par(mar=c(3,7,1,1), xaxs='r', mgp=c(2,1,0))
maxAb <- round(max(log_Ctrl.rel), digits=2)
plot(1, type='n', ylim=c(0.8, (ncol(log_Met.rel)*2)-0.8), xlim=c(0,maxAb), 
     ylab=expression(italic('')), xlab='log Relative Abundance', xaxt='n', yaxt='n', cex.lab=1)
index <- 1

for(i in colnames(log_Met.rel)){
  stripchart(at=index+0.35, log_Met.rel[,i], 
             pch=21, bg='black', method='jitter', jitter=0.15, cex=1, lwd=0.5, add=TRUE)
  stripchart(at=index-0.35, log_Ctrl.rel[,i], 
             pch=21, bg='white', method='jitter', jitter=0.15, cex=1, lwd=0.5, add=TRUE)
  if (i != colnames(log_Met.rel)[length(colnames(log_Met.rel))]){
    abline(h=index+1, lty=2)
  }
  segments(median(log_Met.rel[,i]), index+0.9, median(log_Met.rel[,i]), index+0.1, lwd=2.5, col="black") #adds line for median (can edit to add line for mean instead)
  segments(median(log_Ctrl.rel[,i]), index-0.9, median(log_Ctrl.rel[,i]), index-0.1, lwd=2.5, col="black") #adds line for median (can edit to add line for mean instead)
  index <- index + 2
}
axis(side=1, at=c(0,maxAb), label=c(0,maxAb), cex.axis=1, tck=-0.02)
minors <- c(0.1,0.28,0.44,0.58,0.7,0.8,0.88,0.94,0.98)
axis(side=1, at=minors, label=rep('',length(minors)), tck=-0.01)
legend('bottomright', legend=c('Met PN', 'Ctrl PN'),
       pch=c(21, 21), pt.bg=c('black','white'), bg='white', pt.cex=1, cex=0.8)
####TAXONOMY CODE TO EDIT LABEL NAMES######
###delete extra periods and numbers in taxonomy name####
mini.taxonomy <- as.vector(colnames(log_Met.rel))

#define length of vector for the "for loop"
mini.length <- as.numeric(length(mini.taxonomy))

#loop through names and remove unwanted characters
for(i in 1:mini.length){
  mini.string <- as.character(mini.taxonomy[i])
  mini.string.edit <- removeNumbers(mini.string) #remove numbers
  name <- str_replace_all(mini.string.edit, "[[:punct:]]", "") #removes "(" ")" and all special characters
  mini.taxonomy[i] <- name
}
#apply changes
colnames(log_Met.rel) <- as.character(mini.taxonomy)
######END EDIT LABEL NAMES######
#print labels
axis(2, at=seq(1,index-2,2)+0.6, labels=colnames(log_Met.rel), las=1, line=-0.5, tick=F, cex.axis=0.75, font = 3)
