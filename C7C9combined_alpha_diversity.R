##Code for plotting alpha diversity metrics 
##created by Zach Carlson.  4/13/2018
##code from GitHub repo kassambara/easyGgplot2@master
##from URL https://api.github.com/repos/kassambara/easyGgplot2/zipball/master

##if you don't have easyGgplot2, run the code below##
##CODE:
# install.packages("devtools")
# library(devtools)
# install_github("easyGgplot2", "kassambara")
# library(easyGgplot2)
#####################################################

#load libraries
library(easyGgplot2)
library(ggsignif) #plots significance

# VARIABLES #####################################################################
WORKING.DIRECTORY <- "H:/My Documents/MothurFiles/C7C9_Combined/AlphaDiversity/"
SUMMARY.FILE <- "combined.final.8wk.groups.summary.csv"
DESIGN.FILE <- "combined_8WK.design.txt"
TITLE <- "8WK Pup Shannon Diversity"

#################################################################################

setwd(WORKING.DIRECTORY)
stats <- read.csv(file = SUMMARY.FILE, header = TRUE)
groups <- read.table(file = DESIGN.FILE, header = TRUE, row.names = 1) 

#####Calculate Number of Mets and Ctrls######
groups_ordered<- groups[order(groups$Group, groups$SampleID),] #order table for splittin
all <- split(groups_ordered, groups_ordered$Group) 
met <- all$'MetPN' 
ctrl <- all$'CtrlPN' 
MET.LENGTH <- as.numeric(nrow(met))
CTRL.LENGTH <- as.numeric(nrow(ctrl))
type <-c(rep("Ctrl PN", CTRL.LENGTH), rep("Met PN", MET.LENGTH))
############CALCULATION COMPLETE############

type <- c(rep("Met PN", MET.LENGTH), rep("Ctrl PN", CTRL.LENGTH)) #specify experimental group of each sample

#select shannon data and type and combine
Ratio <- subset(stats, select = shannon)
stats.shannon <- cbind(type, Ratio)

#plot shannon
ggplot2.stripchart(stats.shannon, 
                   xName='type',
                   yName='shannon',
                   size = 10, #size of points
                   addBoxplot = TRUE,
                   mainTitle = TITLE,
                   backgroundColor = "white",
                   fill='#FFAAD4',
                   removePanelGrid = TRUE,
                   removePanelBorder = TRUE,
                   axisLine = c(0.5, "solid", "black"),
                   groupName = 'type',
                   groupColors = c('#999999', '#E69F00'))
#t.test
#t.test(shannon ~ type, data = stats.shannon)

#################################
#plot shannon + significance test#
ggplot(stats.shannon, aes(x=type, y=shannon), mainTitle = TITLE) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Ctrl PN", "Met PN")), 
              map_signif_level=TRUE)
#################################
