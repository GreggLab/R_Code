##Code for a BF ratio. Modeled after code written by Anna Seekatz on her GitHub.com account
##under Rcode_erinsubset/erinsubet_Fig3.heatmap.R and code sent by Kaylie Bullock in Feb 2018

#load required libraries
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(stringr) #loaded to use str_replace_all() which removes all special characters in a string
library(tm) #loaded to use removeNumbers() which removes any number in a string

######ENTER IN VALUES######
TIMEPOINT <- as.character("8WK")
TITLE <- "8WK B:F Ratio"
TIMEPOINT.FILE <- "combined.final.shared"
NAMES.FILE <- "edited.names.axes.csv"
DESIGN.FILE <- "combined_8WK.design.txt"
RATIO.FILE <- "ratio_8wk.csv" #file name for export.  MUST include .csv file extension
TIMEPOINT.SUBSET <- 245:263
###########################
###########################
#TIMEPOINT.SUBSETS
#ratio_table[17:24,]: mom p7
#ratio_table[32:44,]: p17 
#ratio_table[45:65,]: p19
#ratio_table[156:178,]: 4wk
#ratio_table[245:263,]: 8wk
###########################




#set working directory to where your .shared file is and load file
# setwd("~/Downloads/Work/MothurFiles/C7C9_Combined/BF_Ratio") #for mac users
setwd("H://My Documents/MothurFiles/C7C9_Combined/BF_Ratio") #for windows users

#retrieve a subset of shared and save it to new vector "otu"
#"-c" deletes the columns label and numOtus
shared <- read.table(file = TIMEPOINT.FILE, header = TRUE, row.names = 2)
otu <- subset(shared, select =-c(label, numOtus))

#filter out all otus that OVERALL have 200 reads or under. not worth analyzing
otu.filtered <- otu[, which(colSums(otu) >= 200)]  

###############TAXONOMY CODE###############
#add ID names and uses naming mechanism in edited.names.axes.csv
#NAMING MECHANISM:  switched ID names from 251b17 to P17Ctrl1, etc.
#we'll automatically get a sorted table ordered first by time point and then by group.
#edited.names.axes.csv is just a pcoa file where the names were manually changed, the sample order matches.

#use code below, not needed for this run
real.names <- read.csv(file = NAMES.FILE, header = TRUE)
row.names(otu.filtered) <- row.names(real.names) #row names has correct names

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
  phylum <- sapply(strsplit(string,";"), `[`, 2) #should be 2
  pre.pre.name <- str_c(phylum) ###simple version for now
  pre.name <- removeNumbers(pre.pre.name) #remove numbers
  name <- str_replace_all(pre.name, "[[:punct:]]", "") #removes "(" ")" and all special characters
  #name <- removeWords(name.unclassified, "unclassified")  #removes the word unclassified ###doesn't work yet
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
#################END TAXONOMY CODE#################

#data frame now has correctly edited labels and is ordered alphabetically first by timepoint and then by experimental group.  we can now calculate B/F ratios

##Find Bacteriodetes sum

#because we are cycling through the entire data frame, we need a for loop for both rows and columns
#we have a for loop length for columns already, because we used it for the taxonomy changes, but we don't have one for rows
#length of columns is length.data
#define length of rows
length.rows <- as.numeric(dim(otu.filtered)[1])

#create sums column
otu.filtered$BacSum <- rep.int(0, nrow(otu.filtered))
otu.filtered$FirmSum <- rep.int(0, nrow(otu.filtered))

#cycle through array and calculate sums
for(i in 1:length.rows){
  for(j in 1:length.data){
    #if column name is Bacteroidetes and it isn't zero add it to sum
    if(colnames(otu.filtered)[j] == "Bacteroidetes" & otu.filtered[i,j] != 0){
      otu.filtered[i,"BacSum"] <- otu.filtered[i, "BacSum"] + otu.filtered[i,j]
    }
    #if column name is Firmicutes and it isn't zero add it to sum
    if(colnames(otu.filtered)[j] == "Firmicutes" & otu.filtered[i,j] != 0){
      otu.filtered[i,"FirmSum"] <- otu.filtered[i, "FirmSum"] + otu.filtered[i,j]
    }
    else{
    }
  }
}

#calculate ratio
#create ratio column
otu.filtered$ratio <- rep.int(0, nrow(otu.filtered))

for(i in 1:length.rows){
  otu.filtered[i,"ratio"] <- otu.filtered[i,"BacSum"] / otu.filtered[i, "FirmSum"]
}

#pull just the bacteroidetes sum, firmicutes sum, and ratio so it's easier to look at the table
ratio_table <- subset(otu.filtered, select = c("BacSum", "FirmSum", "ratio"))
ratio.point <- as.data.frame(ratio_table[TIMEPOINT.SUBSET,]) #subsets timepoint
write.csv(ratio.point, file = RATIO.FILE) #exports file for excel or prism

#ADD GROUPS TO EACH SAMPLE
groups <- read.table(file = DESIGN.FILE, header = TRUE, row.names = 1) #import group types
# meta <- 
# ratio.point$group <- c(rep("Ctrl PN", 6), rep("Met PN", 9)) #adds group
######

#####Calculate Number of Mets and Ctrls######
groups_ordered<- groups[order(groups$Group, groups$SampleID),] #order table for splittin
all <- split(groups_ordered, groups_ordered$Group) 
met <- all$'MetPN' 
ctrl <- all$'CtrlPN'
MET.LENGTH <- as.numeric(nrow(met))
CTRL.LENGTH <- as.numeric(nrow(ctrl))
ratio.point$group <- c(rep("Ctrl PN", CTRL.LENGTH), rep("Met PN", MET.LENGTH))
############CALCULATION COMPLETE############

#calculate means and se to make plot
ratio.graph <- aggregate(ratio.point$ratio,
                             by = list(ratio.point$group), 
                             FUN = function(x) c(mean = mean(x), sd = sd(x),
                                                 n = length(x)))
ratio.graph <- do.call(data.frame, ratio.graph)
ratio.graph$se <- ratio.graph$x.sd / sqrt(ratio.graph$x.n)
colnames(ratio.graph) <- c("Group", "Ratio", "sd", "n", "se")

plotTop <- 2*max(ratio.graph$Ratio) - (0.5 * min(ratio.graph$Ratio))

#plot
barCenters <- barplot(height = ratio.graph$Ratio,
                      width = 0.25,
                      space = 0.75,
                      xlim = c(0,1),
                      ylim = c(0, plotTop),
                      names.arg = ratio.graph$Group,
                      ylab = "Ratio",
                      col    = c("white", "black"), #makes control bar white and met bar black
                      beside = true, las = 2,
                      cex.names = 0.75, xaxt = "n",
                      main   = TITLE,
                      border = "black", axes = TRUE)

#add error bars
segments(barCenters, ratio.graph$Ratio - ratio.graph$se * 2, barCenters,
         ratio.graph$Ratio + ratio.graph$se * 2, lwd = 1.5)

#add ends of error bars to them
arrows(barCenters, ratio.graph$Ratio - ratio.graph$se * 2, barCenters,
       ratio.graph$Ratio + ratio.graph$se * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

#test for significance
t.test(ratio.point$ratio ~ ratio.point$group)
