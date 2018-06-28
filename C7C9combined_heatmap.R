##Code for a heatmap. Modeled after code written by Anna Seekatz on her GitHub.com account
##under Rcode_erinsubset/erinsubet_Fig3.heatmap.R and code sent by Kaylie Bullock in Feb 2018

#load required libraries
library(RColorBrewer) #gives color schemes
library(ggplot2) #data plotter
library(gplots) #data plotter
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(stringr) #loaded to use str_replace_all() which removes all special characters in a string
library(tm) #loaded to use removeNumbers() which removes any number in a string

#set working directory to where your .shared file is and load file
setwd("~/Downloads/Work/MothurFiles/C7C9_Combined/Heatmap")
shared <- read.table(file = "combined.final.shared", header = TRUE, row.names = 2)


#retrieve a subset of shared and save it to new vector "otu"
#"-c" deletes the columns label and numOtus
otu <- subset(shared, select =-c(label, numOtus))

#filter out all otus that OVERALL have 200 reads or under. not worth analyzing
otu.filtered <- otu[, which(colSums(otu) >= 200)]  

#########################################################################################
#test code for HEATMAP
otu <- subset(shared, select =-c(label, numOtus))
otu.filtered <- subset(otu.filtered, select = c(Otu0008, Otu0003, Otu0020, Otu0041, Otu0017, Otu0073, Otu0033, Otu0018, Otu0037, Otu0089, Otu0032, Otu0025, Otu0048, Otu0091, Otu0058, Otu0031, Otu0072)) #otus selected from LDA values and tests
#########################################################################################


#add ID names and uses naming mechanism in edited.names.axes copy.csv
#NAMING MECHANISM:  switched ID names from 251b17 to P17Ctrl1, etc.
#we'll automatically get a sorted table ordered first by time point and then by group.
#edited.names.axes copy.csv is just another saved version of combined.final.shared where we manually changed all the ID names
real.names <- read.table(file = "edited.names.axes.csv", header = TRUE)
row.names(otu.filtered) <- real.names[,"group"] #second column has the names we are importing

#order table by ID name alphabetically, so we can subset by timepoint later on.  
otu.filtered <- otu.filtered[ order(row.names(otu.filtered)), ]  


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

#export completely edited otu.filtered to view in excel to get row number subsets.  Can view in R but it's a big table and the program might crash!  
#write.csv(otu.filtered, file = "otu.filtered.csv", row.names = TRUE)

#data frame now has correctly edited labels and is ordered alphabetically first by timepoint and then by experimental group.  we can subset based on timepoints and graph just those.
##NOTE:  could simplify with a for loop, etc.  will do eventually..

#separate based on time period
#time periods
#p17: 32:44
#p19: 45:65
#4wk: 156:178
#8wk: 245:263

otu.filtered.complete <- otu.filtered[45:65,] #OTU FOR P19s

#change otu.filtered to a matrix using as.matrix
otu.matrix <- as.matrix(otu.filtered.complete)

#perform relative abundance calculations
otu.rel <- otu.matrix/rowSums(otu.matrix)
otu.rel.max <- apply(otu.rel, 2, max) #2 means columns, apply function max to columns in otu.rel
otu.rel.filtered <- otu.rel[, otu.rel.max>0.02]
row.names(otu.rel.filtered) <- c(rep("Ctrl PN", 13), rep("Met PN", 8))
#create color scheme and relative abundance scales
myCol3 <- colorRampPalette(brewer.pal(8,"GnBu"))(8) #GnBu.. classic green to blue pallette
myBreaks2 <- c(0, 0.001, 0.003, 0.01, 0.05, 0.10, 0.50, 0.80, 1) 

#export tables to analyze for significance in prism!
#write.csv(otu.rel.filtered, file = "8wk.heatmap.csv", row.names = TRUE)

##PROCESSING COMPLETE, PRODUCE YOUR FIGURES!
#heatmap for p19
heatmap.2(otu.rel.filtered, 
          dendrogram = 'none', #removes dendogroms
          key      = TRUE, #removes color key
          #col      = myCol3, #regular heatmap
          #breaks   = myBreaks2,  #regular heatmap
          Colv = TRUE, #doesn't reorder column
          Rowv = FALSE, #doesn't reorder rows so mets and ctrls are separate
          scale    = "column", #used if want red/blue z-distribution
          tracecol = "#303030", #used if want red/blue z-distribution
          col = bluered, #used if want red/blue z-distribution
          main     = "P19 Heatmap", #title
          trace    = 'none',
          cexRow = 1, #sets y axis label text size
          cexCol = 1, #reduces x axis label text size
          margins = c(6,6), #sets graph margins
          srtCol = 30) #adds 30º angle to x axis labels



















