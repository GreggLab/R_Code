#!/usr/bin/env Rscript
#
# This script makes the PCoA plots.
#
# usage:
# 65_make_PCoA_biplot [--no-arrows] <pcoa_axes_file_from_mothur> <loadings_file> [<corr.axes_file> <otu_taxonomy_file>] <design_file> <out_put_filename>
#
# Parameter description:
#
# <pcoa_axes_file_from_mothur>   .pcoa.axes file from mothur   
# <loadings_file>                .pcoa.loadings file rrom mothur
# <corr.axes_file>               .corr.axes file from mothur
# <otu_taxonomy_file>            .cons.taxonomy file from mothur
# <design_file>                  mothur-compatible design file, two columns with header row
# <out_put_filename>             What you want to call the output file (with .pdf suffix)

# whether to draw the arrows (TRUE or FALSE), can be overridden by --no-arrows command line argument
DRAW_ARROWS = TRUE
# size of arrow label text, something between 0.5 and 1.0 or so
CEX_LABEL = 0.65
# LABEL_PADDING
LABEL_PADDING = 1.2
# size of dots, something between 1.0 and 2.0 is reasonable
DOT_SIZE=1.5
# (internal) name of the category when information is missing, needs to correspond to 23_make_design_files
NO_INFO = 'NO_INFO'
# unclassified taxon string
UNCLASSIFIED='unclassified'
# text we will actually print instead of NO_INFO
NO_INFO_TEXT = 'not specified'
# color of the NO_INFO category
NO_INFO_COLOR = 'gray'
# numbers of colors per luminance level
COLS_PER_LUM=5
# hue that the color aligorithm starts with on the color wheel (values 0-360)
START_HUE = 100
# p value cut off for corr.axes arrows
CORR_AXES_P_CUTOFF = 0.05
# max number of arrows drawn
CORR_AXES_MAX_ARROWS = 10

# don't draw arrows of OTUs with low sequence count (FIXME: does this make any sense? If yes, what is a good value here?)
CORR_AXES_MIN_SEQ_COUNT = 5000

library(colorspace, quiet = TRUE)
library(shape)
library(plyr) #tool for splitting, applying, and combining data
library(stringr) #loaded to use str_replace_all() which removes all special characters in a string
library(tm) #loaded to use removeNumbers() which removes any number in a string

#args <- commandArgs(TRUE)
args = vector(mode="character", length=6)
args[1] = 'combined.final.8wk.thetayc.0.03.lt.pcoa.axes'
args[2] = 'combined.final.8wk.thetayc.0.03.lt.pcoa.loadings'
args[3] = 'combined.final.8wk.pearson.corr.axes'
args[4] = 'combined.final.0.03.cons.taxonomy'
args[5] = 'combined_8WK.design.txt'
args[6] = 'PCoA_biplot.pdf'

# setwd("~/Downloads/Work/MothurFiles/C7C9_Combined/PCoA_Biplots") #for mac users
setwd("H:/My Documents/MothurFiles/C7C9_Combined/PCoA_Biplots")#for windows users
##################################################
# (very primitive) command line argument parsing
#
if (args[1] == '--no-arrows') {
  if (length(args) == 5) {
    DRAW_ARROWS = FALSE
    pcoa_axes_file = args[2]
    loadings_file  = args[3]
    design_file    = args[4]
    out_file       = args[5]
  } else {
    cat('expect 4 more arguments in addition to --no-arrows, ', length(args)-1, 'were supplied.\n')
  }
} else {
  if (length(args) == 6) {
    pcoa_axes_file = args[1]
    loadings_file  = args[2]
    corr_axes_file = args[3]
    otu_tax_file   = args[4]
    design_file    = args[5]
    out_file       = args[6]
  } else {
    cat('expect 6 arguments if --no-arrows is not given, ', length(args), 'were supplied.\n')
  }
}#endif

cat('[65_make_PCoA_biplot.R] ')
if (DRAW_ARROWS == FALSE) {
  cat('not drawing arrows ')
}
cat(design_file)

dat = read.table(file = pcoa_axes_file, header=TRUE)
dat = dat[, c('group', 'axis1', 'axis2')  ]
# If group/sample names are made up of digits only, problems occur when these fields are read 
# they'll become integers.
# We need to make sure groups are factors / same goes for the columns of the design files below
dat$group = as.factor(dat$group)
loadings = read.table(file = loadings_file, header=TRUE)
if (DRAW_ARROWS) {
  corr.axes = read.table(file = corr_axes_file, header=TRUE, row.names = 1)
  otu_taxonomy = read.table(file = otu_tax_file, header = TRUE, row.names = 1)
}
design = read.table(file = design_file, col.names=c('group', 'category'), header = TRUE, colClasses = c('factor'))
# make sure category levels are ordered as they come in from the design file,
# this is the order they will be printed in the legend
# and avoids the default alphabetical? level order
design$category <- factor(design$category, levels = unique(design$category))

if (DRAW_ARROWS) {
  ######################################################
  #
  # Prepare corr.axes data
  #
  # add taxonomic classification
  corr.axes = merge(corr.axes, otu_taxonomy[rownames(corr.axes),], by='row.names', all=TRUE)
  # fix rownames (as merge introduces a "Row.names" column (as first column) instead of proper row.names)
  rownames(corr.axes) = corr.axes$Row.names
  corr.axes = corr.axes[,2:8] #remove redundant column
  # pick lowest-ranked classified taxon from taxonomy column
  # assume its formatted like: Bacteria(100);some_phylumn(100);"blabla"(99);unclassified(100);unclassified(100);unclassified(100);
  taxonomy <- as.vector(corr.axes[,"Taxonomy"])
  
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
  corr.axes$Taxonomy <- taxonomy
  # sort by length, pick first 5, used in biplot later
  corr.axes = corr.axes[ order(-corr.axes$length), ]
  # remove arrows with high p values for both axes
  corr.axes = corr.axes[ corr.axes$p.value < CORR_AXES_P_CUTOFF | corr.axes$p.value.1 < CORR_AXES_P_CUTOFF  ,]
  # remove arrows based on low sequence count
  corr.axes = corr.axes[ corr.axes$Size > CORR_AXES_MIN_SEQ_COUNT,]
  # Limit number of arrows drawn
  if (dim(corr.axes)[1] > CORR_AXES_MAX_ARROWS) {
    corr.axes = corr.axes[1:CORR_AXES_MAX_ARROWS,]
  }
}#endif

# Get the coordinate range for the plotting window to get a design independend common range
# as we may lose that info after the merge() call
xlim = c(min(dat$axis1), max(dat$axis1))
ylim = c(min(dat$axis2), max(dat$axis2))

# assign category to each datapoint:
# merge order, sort: we want to keep the order of design
dat = merge(design, dat, by='group', sort = FALSE)

###########################################chroma
#
# SHAPES
pch <- as.vector(c(16, 1)) #shape, 16 = filled dot, 1 = open dot
#finally, assign the shapes to the data
data_pch <- as.vector(pch[dat$category])

################################################3
#
# PLOT
#

#prepare titles
info=design_file
info = sub('.design.txt','', info)
info=strsplit(info, '___')[[1]]
if (length(info) == 1) {
  # copy main title to legend title for design file names not conforming to the ___ syntax
  info = append(info, info[1])
}
#main = info[1]
#if (main == '') {
#    main = 'all samples'
#}
#main = gsub('__', ', ', main)
#main = gsub('_', ' ', main)
#sub = paste('AMOVA test results are in file [...].', design_file, '.amova', sep='')
# graphics parameters
#pdf(file = out_file,  11, 8.5, onefile = TRUE)
par(mai=c(1,0.9,0.3,1.9), cex.sub=0.75)
par(xpd = TRUE, col.axis = 'gray77', fg='gray77', bty='n')
#par(xpd=TRUE) # allow legend to be outside the plot
x_label = paste('Axis 1 (', signif(loadings$loading[1], digits=3), '%)', sep='')
y_label = paste('Axis 2 (', signif(loadings$loading[2], digits=3), '%)', sep='')

# set up everything, but don't draw any dots yet (because arrows come first),
plot(dat$axis1, dat$axis2, type='n', xlab=x_label, ylab=y_label, cex=DOT_SIZE, xlim=xlim, ylim=ylim)

if (DRAW_ARROWS) {
  ###########################
  # Draw corr.axes arrows
  #
  # calculate scaling factor for arrows (so that they fit into the plot)
  scale_north = max(dat$axis2) / max(corr.axes$axis2)
  scale_east  = max(dat$axis1) / max(corr.axes$axis1)
  scale_south = min(dat$axis2) / min(corr.axes$axis2)
  scale_west  = min(dat$axis1) / min(corr.axes$axis1)
  scale = min(abs(c(scale_north, scale_east, scale_south, scale_west)))
  # draw the arrows
  if (DRAW_ARROWS) {
    Arrows(0, 0, corr.axes$axis1*scale, corr.axes$axis2*scale, col='gray')
  }
}

# draw the dots
points(dat$axis1, dat$axis2, col="black", pch=data_pch, cex=DOT_SIZE)

################################
# prepare for legend:
#
par(fg='black')
# get categories for legend
categories = as.vector(c("Met PN", "Ctrl PN"))
title = info[2]
title = gsub('__',', ', title)
title = gsub('_',' ', title)
# capitalize
#title = paste(toupper(substring(title,1,1)), substring(title, 2), sep='')
#
# legend positioning !@$%@!#$:
#   mai specifies margins to make place on the right side
#   xpd=TRUE (or was it NA?) to allow the legend to shine though outside the plot
#
#   We calculate the legend position in term of plot data point, as that seems to
#   the way to go to get some consitent position for any data.  Just using the maximum
#   data values for x and y will put the legend in the top-right corner but we need
#   to move it a bit further to the right, about 1/2 inch on the letter sized PDF:
plot_width_in = 7.5 #estimated width of plotted region in inch
adj_in = 0.5 # legend x-adjustment in inch
plot_width_data = xlim[2] - xlim[1]
adj_x = plot_width_data * adj_in / plot_width_in
# finally set legend position
x = as.vector(xlim[2] + adj_x)
y = as.vector(ylim[2])

legend(x = x, y = y,
       legend=categories,
       #inset=c(-0.24,0),
       pch=pch,
       col="black",
       bty="n",
       cex=0.8,
       title = title,
       title.adj = 0,
       text.width=0.25,
       xpd=TRUE)


if (DRAW_ARROWS) {
  ###########################
  # Arrow labels
  #
  # arrow label scaling factor to print test a bit away from arrow tip
  label_offset = (xlim[2] - xlim[1])*0.05
  slopes = abs(corr.axes$axis2 / corr.axes$axis1)
  x_offsets = label_offset / sqrt(slopes^2+1) * sign(corr.axes$axis1)
  y_offsets = slopes * x_offsets * sign(corr.axes$axis1) * sign(corr.axes$axis2)
  # fix y_offset for arrows pointing straight up or down (infinite slope!) y_offset: NaN --> label_offset
  # sign depends of sign of Inf slope
  y_offsets[is.nan(y_offsets)] = sign(slopes[is.infinite(slopes)]) * label_offset
  ######################################
  #calculate label positions
  #
  # replace arrow coordinates with label coordinates
  corr.axes$axis1 = corr.axes$axis1 * scale + x_offsets
  corr.axes$axis2 = corr.axes$axis2 * scale + y_offsets
  # build label from Otu and taxon and add a label column
  rownames(corr.axes) = sub('Otu[0]*', 'OTU ', rownames(corr.axes))
  label = paste(rownames(corr.axes), corr.axes$Taxonomy, sep='\n')
  corr.axes = cbind(corr.axes, label)
  #########################################################
  # shift labels around to eliminate overlapping
  #
  # split data at y=0 and order so we pull overlapping labels away from the plot origin (up or down)
  # labels won't be moved left-right since they are much wider than tall
  upper_half = corr.axes[corr.axes$axis2 >= 0,]
  upper_half = upper_half[ order(-upper_half$axis2),]
  lower_half = corr.axes[corr.axes$axis2 <  0,]
  lower_half = lower_half[ order(lower_half$axis2),]
  overlap = TRUE
  # the while loop *should* terminate on its own, well before MAX_ITER loops for low number of arrows, MAX_ITER is just a safeguard
  MAX_ITER = 50
  iter = 0
  # run iteration until no overlap
  while (overlap && iter < MAX_ITER) {
    iter = iter + 1
    overlap = FALSE
    
    # upper half
    for (i in seq(length.out = max(0, length(rownames(upper_half))-1))) {
      # if upper half has more than one row, enter this loop for every row
      # except the last one
      # i is the index of the label we want to adjust (if necessary)
      for ( j in seq(i+1, length(rownames(upper_half))) ) {
        # assume upper_half has at least two rows
        # j is the index of the label we compare label i with
        # j runs from i+1 to last index of upper half
        # calculate distance:
        x_dist = upper_half[i,'axis1'] - upper_half[j,'axis1']
        y_dist = upper_half[i,'axis2'] - upper_half[j,'axis2']
        # width and height of this and the next label
        width       = strwidth( upper_half[i,'label'], cex = CEX_LABEL * LABEL_PADDING)
        next_width  = strwidth( upper_half[j,'label'], cex = CEX_LABEL * LABEL_PADDING)
        height      = strheight(upper_half[i,'label'], cex = CEX_LABEL * LABEL_PADDING)
        next_height = strheight(upper_half[j,'label'], cex = CEX_LABEL * LABEL_PADDING)
        # do they overlap?
        x_overlap = abs(x_dist) < (width  + next_width ) / 2
        y_overlap = abs(y_dist) < (height + next_height) / 2
        ### BEGIN DEBUG CODE ###
        if (is.na(x_overlap) || is.na(y_overlap)) {
          print('!error condition! -- upper half de-overlap')
          print(upper_half)
          cat('length: ', length(upper_half), '\n')
        }
        #### END DEBUG CODE ####
        if (x_overlap && y_overlap) {
          # if overlap calculate and apply adjustement 
          label_adj = (height + next_height) / 2 - y_dist
          # * 1.00001 avoid possible infinite loop
          upper_half[i, 'axis2'] = upper_half[i, 'axis2'] + label_adj * 1.00001
          #cat('Adjusting labels', rownames(upper_half)[i], '<--', rownames(upper_half)[j], label_adj, '\n')
          # assume there's more overlap (or we just introduced some)
          overlap = TRUE
          break
        }
      }#end for
    }#end for
    
    # lower half
    for (i in seq(length.out = max(0, length(rownames(lower_half))-1))) {
      # if lower half has more than one row, enter this loop for every row
      # except the last one
      # i is the index of the label we want to adjust (if necessary)
      for ( j in seq(i+1, length(rownames(lower_half))) ) {
        # assume lower_half has at least two rows
        # j is the index of the label we compare label i with
        # j runs from i+1 to last index of lower half
        # calculate distance, y-dist will be negative:
        x_dist = lower_half[i,'axis1'] - lower_half[j,'axis1']
        y_dist = lower_half[i,'axis2'] - lower_half[j,'axis2']
        # width and height of this and the next label
        width       = strwidth( lower_half[i, 'label'], cex = CEX_LABEL * LABEL_PADDING)
        next_width  = strwidth( lower_half[j, 'label'], cex = CEX_LABEL * LABEL_PADDING)
        height      = strheight(lower_half[i, 'label'], cex = CEX_LABEL * LABEL_PADDING)
        next_height = strheight(lower_half[j, 'label'], cex = CEX_LABEL * LABEL_PADDING)
        # do they overlap?
        x_overlap = abs(x_dist) < (width  + next_width ) / 2
        y_overlap = abs(y_dist) < (height + next_height) / 2
        ### BEGIN DEBUG CODE ###
        if (is.na(x_overlap) || is.na(y_overlap)) {
          print('!error condition! -- lower half de-overlap')
          print(lower_half)
          cat('length: ', length(lower_half), '\n')
        }
        #### END DEBUG CODE ####
        if (x_overlap && y_overlap) {
          # if overlap calculate and apply adjustement 
          label_adj = (height + next_height) / 2 + y_dist
          # * 1.00001 avoid possible infinite loop
          lower_half[i, 'axis2'] = lower_half[i, 'axis2'] - label_adj * 1.00001
          #cat('Adjusting labels', rownames(lower_half)[i], '<--', rownames(lower_half)[j], -label_adj, '\n')
          # assume there's more overlap (or we just introduced some)
          overlap = TRUE
          break
        }
      }#end for
    }#end for
  }# end while
  #cat('Needed', iter, 'label adjustment iterations to avoid overlapping.\n')
  
  # put halfs back together
  corr.axes = rbind(upper_half, lower_half)
  
  # print the labels
  par(lheight = 0.8)
  text(corr.axes$axis1,
       corr.axes$axis2,
       label   = corr.axes$label,
       cex     = CEX_LABEL
  )
}#endif DRAW_ARROWS

cat(' done.\n')


