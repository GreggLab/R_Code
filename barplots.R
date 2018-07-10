#Stacked Community Plots
#6-19-18
#Adapted from Christine Bassis, Ph. D. code
#Adapted by Zach Carlson, email: zcarlson@umich.edu
library(RColorBrewer) #gives color schemes
library(ggplot2) #data plotter
library(gplots) #data plotter
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(stringr) #loaded to use str_replace_all() which removes all special characters in a string
library(tm) #loaded to use removeNumbers() which removes any number in a string
library(shiny)

######ENTER IN VALUES######
TIMEPOINT <- as.character("8WK")
TITLE <- "8WK Relative Abundances"
TAXONOMY.FILE <- "combined.final.0.03.cons.taxonomy"
TIMEPOINT.FILE <- "combined.final.8wk.shared"
DESIGN.FILE <- 'combined_8WK.design.txt'#Cols: Row, SampleID, Group
PHYLA <- c("Actinobacteria", "Bacteria", "Bacteroidetes", "CandidatusSaccharibacteria","Chloroflexi", "Deferribacteres", "Firmicutes", "Proteobacteria", "Tenericutes", "Verrucomicrobia")
REL.FILE <- "otu.rel.t."
TXT <- ".txt"
TAXONOMY.COLOR.FILE <- "taxonomy.color."
OTU.COMPLETE <- "otu.complete."
###########################
#get into correct directory
# setwd("~/Downloads/Work/MothurFiles/C7C9_Combined/Stacked_Community_Plots") #mac users
setwd("H:/My Documents/MothurFiles/C7C9_Combined/Stacked_Community_Plots") #for windows users

#read in raw files 
tax <- read.table(file=TAXONOMY.FILE, 
                       row.names = 1,
                       header=TRUE,
                       check.names=FALSE,
                       comment.char="") #will become otu.class
otu.info <- read.table(file=TIMEPOINT.FILE, header=TRUE, row.names = 2) #will become rm_g
meta <- read.table(file=DESIGN.FILE, row.names = 1, header =TRUE)

#get length values to use later on
SAMPLE.LENGTH <- as.numeric(nrow(meta))
TAX.SAMPLE.LENGTH <- SAMPLE.LENGTH + 4 #adds OTU, phyla, genus, color
TAXONOMY.LENGTH <- as.numeric(nrow(tax))
otu.info <- subset(otu.info, select = -c(label, numOtus))
PHYLA.LENGTH <- length(PHYLA)

###CREATE OTU.CLASS FILE###
# - pull phlyum
# - pull genus
##following code changes label "OTU0001" to corresponding genus/phyla.##

#pull just taxonomy string 
taxonomy <- as.vector(tax[,"Taxonomy"])

#create table to store phyla and genus values
otu.class <- data.frame(matrix(NA, nrow = TAXONOMY.LENGTH, ncol = 5)) #5 for OTU, Size, Phylum, Genus, Full
colnames(otu.class) <- c("OTU", "Size", "Phylum", "Genus", "Full")
otu.class$OTU <- rownames(tax) #add otu numbers
otu.class$Full <- tax$Taxonomy #add taxonomy
otu.class$Size <- tax$Size #add size

#changes lengthy taxonomic info to just a string of the phlya and genus
#sapply(): applys a function across a string
#strsplit(): splits string based on the character ";" and selects the 3rd value or 4th value (phlya and genus respectively)
for(i in 1:TAXONOMY.LENGTH){
  string <- as.character(taxonomy[i]) #maintains character status
  phyla <- sapply(strsplit(string,";"), `[`, 2) #should be 2
  genus <- sapply(strsplit(string,";"), `[`, 6) #should be 6
  phyla <- removeNumbers(phyla) #remove numbers
  genus <- removeNumbers(genus) #remove numbers
  phyla <- str_replace_all(phyla, "[[:punct:]]", "") #removes "(" ")" and all special characters
  genus <- str_replace_all(genus, "[[:punct:]]", "") #removes "(" ")" and all special characters
  phyla <- str_replace_all(phyla, "unclassified", "") #removes "unclassified"
  genus <- str_replace_all(genus, "unclassified", "") #removes "unclassified"
  otu.class[i,"Phylum"] <- phyla
  otu.class[i,"Genus"] <- genus
}
###OTU.CLASS IS COMPLETED###

###CREATE RM_G FILE###
###PART 1. CREATE META COLOR FILE###
taxonomy.color <- data.frame(matrix(NA, nrow = TAXONOMY.LENGTH, ncol = 4)) #4 for OTU, Phylum, Genus, Color
colnames(taxonomy.color) <- c("OTU", "Phylum", "Genus", "Color")
taxonomy.color$OTU <- rownames(tax) #add otu numbers
taxonomy.color$Phylum <- otu.class$Phylum
taxonomy.color$Genus <- otu.class$Genus
colors <- as.vector(colorRampPalette(brewer.pal(9,"Pastel1"))(9))
colors[10] <- '#ef8f8d' #add a tenth color becaue pastel1 only offers 9

#ran table(taxonomy.color$Phylum) to get all different phylum present
#and in what frequency
#created color palette with as many colors as different phylum present
#assign colors based on phylum, non-conforming samples are assigned black
for(i in 1:TAXONOMY.LENGTH){
  if(taxonomy.color$Phylum[i] == 'Actinobacteria'){
    taxonomy.color$Color[i] <- colors[10]
  }
  else if(taxonomy.color$Phylum[i] == 'Bacteria'){
    taxonomy.color$Color[i] <- colors[9]
  }
  else if(taxonomy.color$Phylum[i] == 'Bacteroidetes'){
    taxonomy.color$Color[i] <- colors[8]
  }
  else if(taxonomy.color$Phylum[i] == 'CandidatusSaccharibacteria'){
    taxonomy.color$Color[i] <- colors[7]
  }
  else if(taxonomy.color$Phylum[i] == 'Chloroflexi'){
    taxonomy.color$Color[i] <- colors[6]
  }
  else if(taxonomy.color$Phylum[i] == 'Deferribacteres'){
    taxonomy.color$Color[i] <- colors[5]
  }
  else if(taxonomy.color$Phylum[i] == 'Firmicutes'){
    taxonomy.color$Color[i] <- colors[4]
  }
  else if(taxonomy.color$Phylum[i] == 'Proteobacteria'){
    taxonomy.color$Color[i] <- colors[3]
  }
  else if(taxonomy.color$Phylum[i] == 'Tenericutes'){
    taxonomy.color$Color[i] <- colors[2]
  }
  else if(taxonomy.color$Phylum[i] == 'Verrucomicrobia'){
    taxonomy.color$Color[i] <- colors[1]
  }
  else{
    taxonomy.color$Color[i] <- '#000000' #assigns black to unclassified
  }
}
###PART 1. META COLOR FILE IS COMPLETE###

###PART 2. CREATE RELATIVE ABUNDANCE FILE###
otu.matrix <- as.matrix(otu.info)
otu.rel <- otu.matrix/rowSums(otu.matrix) #check with 'rowSums(otu.rel)' all rows = 1
#otu.rel <- subset(otu.rel, select = -c(label, numOtus))
otu.rel.t <- t(otu.rel) #transpose #check with 'colSums(otu.rel.t)' all columns = 1
otu.rel.t <- as.data.frame(otu.rel.t) #make it so we can add back in OTU names without changing to list
otu.rel.t$OTU <- rownames(otu.rel.t) #add OTUs column so we can merge
rm_g <- merge(taxonomy.color, otu.rel.t, by.x = "OTU", by.y = "OTU")
###PART 2. COMPLETE###
### rm_g FILE COMPLETED ###



otubar <- as.matrix(subset(rm_g, select =-c(Genus, Color, OTU))) #delete all columns
rownames(otubar) <- otubar[,"Phylum"] #keep phylum names by saving them as rownames
otubar <- subset(otubar, select = -c(Phylum)) #delete phylum column so all cells are numeric
barg <- as.data.frame(t(otubar))
barg$SampleID <- meta$SampleID #add IDs
barg$Group <- meta$Group #add groups
col.gen <- as.character(rm_g$Color)
bar_ordered<- barg[order(barg$Group, barg$SampleID),] #order table for splitting

#splits mets and controls
all <- split(bar_ordered, bar_ordered$Group) 
met <- all$'MetPN' 
ctrl <- all$'CtrlPN'
MET.LENGTH <- as.numeric(nrow(met))
CTRL.LENGTH <- as.numeric(nrow(ctrl))
FINAL.LENGTH <- as.numeric(ncol(met))

###MAKE MET FILE###
barmet <- subset(met, select = -c(SampleID, Group))
barmet <- as.matrix(barmet)
class(barmet) <- "numeric" #change matrix to numeric form rather than character
colnames(barmet) <- rm_g[,"Phylum"] #removes numbers from colnames. e.g. Bacteroidetes.1 -> Bacteroidetes

###phyla only distribution###
rows.length <- as.numeric(nrow(barmet))
cols.length <- as.numeric(ncol(barmet))

#make sum columns
barmet <- as.data.frame(barmet)
barmet$ActinobacteriaSums <- rep.int(0, nrow(barmet))
barmet$BacteriaSums <- rep.int(0, nrow(barmet))
barmet$BacteroidetesSums <- rep.int(0, nrow(barmet))
barmet$CandidatusSaccharibacteriaSums <- rep.int(0, nrow(barmet))
barmet$ChloroflexiSums <- rep.int(0, nrow(barmet))
barmet$DeferribacteresSums <- rep.int(0, nrow(barmet))
barmet$FirmicutesSums <- rep.int(0, nrow(barmet))
barmet$ProteobacteriaSums <- rep.int(0, nrow(barmet))
barmet$TenericutesSums <- rep.int(0, nrow(barmet))
barmet$VerrucomicrobiaSums <- rep.int(0, nrow(barmet))

for(i in 1:rows.length){
  for(j in 1:cols.length){
    if(colnames(barmet)[j] == 'Actinobacteria'){
      barmet$ActinobacteriaSums[i] <- barmet[i,j] + barmet$ActinobacteriaSums[i]
    }
    else if(colnames(barmet)[j] == 'Bacteria'){
      barmet$BacteriaSums[i] <- barmet[i,j] + barmet$BacteriaSums[i]
    }
    else if(colnames(barmet)[j] == 'Bacteroidetes'){
      barmet$BacteroidetesSums[i] <- barmet[i,j] + barmet$BacteroidetesSums[i]
    }
    else if(colnames(barmet)[j] == 'CandidatusSaccharibacteria'){
      barmet$CandidatusSaccharibacteriaSums[i] <- barmet[i,j] + barmet$CandidatusSaccharibacteriaSums[i]
    }
    else if(colnames(barmet)[j] == 'Chloroflexi'){
      barmet$ChloroflexiSums[i] <- barmet[i,j] + barmet$ChloroflexiSums[i]
    }
    else if(colnames(barmet)[j] == 'Deferribacteres'){
      barmet$DeferribacteresSums[i] <- barmet[i,j] + barmet$DeferribacteresSums[i]
    }
    else if(colnames(barmet)[j] == 'Firmicutes'){
      barmet$FirmicutesSums[i] <- barmet[i,j] + barmet$FirmicutesSums[i]
    }
    else if(colnames(barmet)[j] == 'Proteobacteria'){
      barmet$ProteobacteriaSums[i] <- barmet[i,j] + barmet$ProteobacteriaSums[i]
    }
    else if(colnames(barmet)[j] == 'Tenericutes'){
      barmet$TenericutesSums[i] <- barmet[i,j] + barmet$TenericutesSums[i]
    }
    else if(colnames(barmet)[j] == 'Verrucomicrobia'){
      barmet$VerrucomicrobiaSums[i] <- barmet[i,j] + barmet$VerrucomicrobiaSums[i]
    }
  }
}

#add sums to new table
sums <- data.frame(matrix(NA, nrow = MET.LENGTH, ncol = PHYLA.LENGTH))
colnames(sums) <- PHYLA
rownames(sums) <- rownames(barmet)
sums$Actinobacteria <- barmet$ActinobacteriaSums
sums$Bacteria <- barmet$BacteriaSums
sums$Bacteroidetes <- barmet$BacteroidetesSums
sums$CandidatusSaccharibacteria <- barmet$CandidatusSaccharibacteriaSums 
sums$Chloroflexi <- barmet$ChloroflexiSums
sums$Deferribacteres <- barmet$DeferribacteresSums
sums$Firmicutes <- barmet$FirmicutesSums
sums$Proteobacteria <- barmet$ProteobacteriaSums
sums$Tenericutes <-barmet$TenericutesSums 
sums$Verrucomicrobia <- barmet$VerrucomicrobiaSums

#transpose for graphing
sums <- as.data.frame(sums)
sums.t <- t(sums)
class(sums.t) <- "numeric"
###FINISH MAKING MET FILE###

###MAKE CTRL FILE###
barctrl <- subset(ctrl, select = -c(SampleID, Group))
barctrl <- as.matrix(barctrl)
class(barctrl) <- "numeric" #change matrix to numeric form rather than character
colnames(barctrl) <- rm_g[,"Phylum"] #removes numbers from colnames. e.g. Bacteroidetes.1 -> Bacteroidetes

###phyla only distribution###
rows.length <- as.numeric(nrow(barctrl))
cols.length <- as.numeric(ncol(barctrl))

#make sum columns
barctrl <- as.data.frame(barctrl)
barctrl$ActinobacteriaSums <- rep.int(0, nrow(barctrl))
barctrl$BacteriaSums <- rep.int(0, nrow(barctrl))
barctrl$BacteroidetesSums <- rep.int(0, nrow(barctrl))
barctrl$CandidatusSaccharibacteriaSums <- rep.int(0, nrow(barctrl))
barctrl$ChloroflexiSums <- rep.int(0, nrow(barctrl))
barctrl$DeferribacteresSums <- rep.int(0, nrow(barctrl))
barctrl$FirmicutesSums <- rep.int(0, nrow(barctrl))
barctrl$ProteobacteriaSums <- rep.int(0, nrow(barctrl))
barctrl$TenericutesSums <- rep.int(0, nrow(barctrl))
barctrl$VerrucomicrobiaSums <- rep.int(0, nrow(barctrl))

for(i in 1:rows.length){
  for(j in 1:cols.length){
    if(colnames(barctrl)[j] == 'Actinobacteria'){
      barctrl$ActinobacteriaSums[i] <- barctrl[i,j] + barctrl$ActinobacteriaSums[i]
    }
    else if(colnames(barctrl)[j] == 'Bacteria'){
      barctrl$BacteriaSums[i] <- barctrl[i,j] + barctrl$BacteriaSums[i]
    }
    else if(colnames(barctrl)[j] == 'Bacteroidetes'){
      barctrl$BacteroidetesSums[i] <- barctrl[i,j] + barctrl$BacteroidetesSums[i]
    }
    else if(colnames(barctrl)[j] == 'CandidatusSaccharibacteria'){
      barctrl$CandidatusSaccharibacteriaSums[i] <- barctrl[i,j] + barctrl$CandidatusSaccharibacteriaSums[i]
    }
    else if(colnames(barctrl)[j] == 'Chloroflexi'){
      barctrl$ChloroflexiSums[i] <- barctrl[i,j] + barctrl$ChloroflexiSums[i]
    }
    else if(colnames(barctrl)[j] == 'Deferribacteres'){
      barctrl$DeferribacteresSums[i] <- barctrl[i,j] + barctrl$DeferribacteresSums[i]
    }
    else if(colnames(barctrl)[j] == 'Firmicutes'){
      barctrl$FirmicutesSums[i] <- barctrl[i,j] + barctrl$FirmicutesSums[i]
    }
    else if(colnames(barctrl)[j] == 'Proteobacteria'){
      barctrl$ProteobacteriaSums[i] <- barctrl[i,j] + barctrl$ProteobacteriaSums[i]
    }
    else if(colnames(barctrl)[j] == 'Tenericutes'){
      barctrl$TenericutesSums[i] <- barctrl[i,j] + barctrl$TenericutesSums[i]
    }
    else if(colnames(barctrl)[j] == 'Verrucomicrobia'){
      barctrl$VerrucomicrobiaSums[i] <- barctrl[i,j] + barctrl$VerrucomicrobiaSums[i]
    }
  }
}

sums.ctrl <- data.frame(matrix(NA, nrow = CTRL.LENGTH, ncol = PHYLA.LENGTH))
colnames(sums.ctrl) <- PHYLA
rownames(sums.ctrl) <- rownames(barctrl)
sums.ctrl$Actinobacteria <- barctrl$ActinobacteriaSums
sums.ctrl$Bacteria <- barctrl$BacteriaSums
sums.ctrl$Bacteroidetes <- barctrl$BacteroidetesSums
sums.ctrl$CandidatusSaccharibacteria <- barctrl$CandidatusSaccharibacteriaSums 
sums.ctrl$Chloroflexi <- barctrl$ChloroflexiSums
sums.ctrl$Deferribacteres <- barctrl$DeferribacteresSums
sums.ctrl$Firmicutes <- barctrl$FirmicutesSums
sums.ctrl$Proteobacteria <- barctrl$ProteobacteriaSums
sums.ctrl$Tenericutes <-barctrl$TenericutesSums 
sums.ctrl$Verrucomicrobia <- barctrl$VerrucomicrobiaSums

sums.ctrl <- as.data.frame(sums.ctrl)
sums.ctrl.t <- t(sums.ctrl)
class(sums.ctrl.t) <- "numeric"

sums.total <- cbind(sums.t, sums.ctrl.t)
colnames(sums.total) <- c(rep("Met PN", MET.LENGTH), rep("Ctrl PN", CTRL.LENGTH))
# graphing both sets:
par(mfrow=c(1,1)) #creates a plot of two columns and one row
par(mar=c(3.3,3,2,1))
par(xpd=T)

barplot(sums.total, 
        las=2, 
        main=TITLE, 
        ylab="Relative Abundance", 
        cex.names=.8, 
        ylim=c(0,1), 
        col=colors, 
        xlim=c(0,40),
        space=c(rep(0.2, MET.LENGTH), 1.5, rep(0.2, CTRL.LENGTH-1)))
legend.x <- MET.LENGTH + CTRL.LENGTH + 5
legend(legend.x, 1,
       legend=rownames(sums.t),
       col=colors,
       fill=colors,
       cex=1, 
       bty="n", 
       ncol=1)		
# 
# not.preg_g<-as.matrix(t(barnot.preg))
# par(mar=c(5,4,2,5))
# par(xpd=T)
# barplot(not.preg_g, las=2, main="Not Pregnant", ylab="Relative abundance", cex.names=0.8, ylim=c(0,1), col=col.gen, xlim=c(0,100))
# ###Note: Had to mess with legend position (check position with locator before plotting 2nd bar plot), when exported all colors/species didn't show up, ended up combining from a couple tries
