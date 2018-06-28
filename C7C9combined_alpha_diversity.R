##Code for plotting alpha diversity metrics 
##created by Zach Carlson.  4/13/2018
##code from GitHub repo kassambara/easyGgplot2@master
##from URL https://api.github.com/repos/kassambara/easyGgplot2/zipball/master
##if you don't have easyGgplot2, run the code below
##CODE:
# install.packages("devtools")
# library(devtools)
# install_github("easyGgplot2", "kassambara")
# library(eastGgplot2)

library(easyGgplot2)
library(ggsignif) #plots significance
setwd("~/Downloads/Work/MothurFiles/C7C9_Combined/AlphaDiversity")
stats <- read.csv(file = "combined.final.p19.groups.summary.csv", header = TRUE)

#specify experimental group of each sample
type <- c("Met PN", "Met PN", "Met PN", "Met PN", "Met PN", "Met PN", "Met PN", "Met PN", "Ctrl PN", "Ctrl PN", "Ctrl PN", "Ctrl PN", "Ctrl PN", "Ctrl PN", "Ctrl PN", "Ctrl PN", "Ctrl PN", "Ctrl PN", "Ctrl PN", "Ctrl PN", "Ctrl PN")

#select shannon data and type and combine
Ratio <- subset(stats, select = shannon)
stats.shannon <- cbind(type, shannon)

#plot shannon
ggplot2.stripchart(stats.shannon, 
                   xName='type',
                   yName='shannon',
                   size = 10, #size of points
                   addBoxplot = TRUE,
                   mainTitle = "P19 Pup Shannon Diversity",
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
ggplot(stats.shannon, aes(x=type, y=shannon), mainTitle = "P17 Shannon Diversity") + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Ctrl", "Met")), 
              map_signif_level=TRUE)
#################################

#select inv simpson data alongside type data and combine
invsimpson <- subset(stats, select = invsimpson)
stats.invsimpson <- cbind(type, invsimpson)

#plot inverse simpson
ggplot2.stripchart(stats.invsimpson, 
                   xName='type',
                   yName='invsimpson',
                   size = 3, #size of points
                   addBoxplot = TRUE,
                   mainTitle = "Mom P7 C7/9PN Inv Simpson Diversity",
                   backgroundColor = "white",
                   fill='#FFAAD4',
                   removePanelGrid = TRUE,
                   removePanelBorder = TRUE,
                   axisLine = c(0.5, "solid", "black"),
                   groupName = 'type',
                   groupColors = c('#999999', '#E69F00'))

#################################
#plot shannon + significance test#
ggplot(stats.invsimpson, aes(x=type, y=invsimpson), mainTitle = "Inv Simpson Diversity") + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Ctrl", "Met")), 
              map_signif_level=TRUE)
#################################