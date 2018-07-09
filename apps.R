library(shiny) 
library(dplyr)
library(DT) #renderDataTable()
library(RColorBrewer) #gives color schemes
library(ggplot2) #data plotter
library(gplots) #data plotter
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(stringr) #loaded to use str_replace_all() which removes all special characters in a string
library(tm) #loaded to use removeNumbers() which removes any number in a string
library(dplyr)
library(readr)
library(rsconnect)
setwd("~/Downloads/Work/MothurFiles/C7C9_Combined/Stacked_Community_Plots")
# rsconnect::setAccountInfo(name='gregglab', token='720B02EC7BE6C0476E21919DBC205C2D', secret='dZmyXnLzIElZ6xctWR9r55ttw73Vu0n2fN1fqvVo')
# rsconnect::deployApp
ui <- fluidPage(
  
  #Application title
  titlePanel("Microbiome Analysis:  Stacked Community Barplots"),
  
  sidebarLayout(

    #Sidebar with input options
    sidebarPanel(
      fileInput(inputId = "taxonomy_file", 
                label = h3("Select Taxonomy File"),
                placeholder = "combined.final.0.03.cons.taxonomy",
                accept = c(".taxonomy", ".csv")),
      checkboxInput("show_taxonomy_table", "Display Taxonomy Table (Allow 1 minute to appear)", FALSE),
      tags$hr(),
      fileInput(inputId = "samples_file", 
                label = h3("Select Shared File"),
                placeholder =  "combined.final.8wk.shared",
                accept = c(".shared", ".csv")),
      checkboxInput("show_samples_table", "Display Samples Table (Allow 1 minute to appear)", FALSE),
      tags$hr(),
      fileInput(inputId = "design_file", 
                label = h3("Select Design File"),
                placeholder = "combined_8WK.design.txt",
                accept = "design.txt"),
      checkboxInput("show_design_table", "Display Design Table (Allow 1 minute to appear)", FALSE),
      tags$hr(),
      textInput(inputId = "timepoint",
                label = h2("Enter in Timepoint"),
                value = "8WK"),
      textInput(inputId = "title",
                label = h2("Enter in Title"),
                value = "8WK Relative Abundances"),
      dateInput("date", 
                h3("Enter Date"),
                format = "mm-dd-yyyy",
                weekstart = 0,
                autoclose = TRUE),
      radioButtons(inputId = "filetype",
                   label = "Select filetype:",
                   choices = c("csv", "tsv"),
                   selected = "csv"),
      selectInput("color_list",
                  label = h2("Select a Color Palette:"),
                  selected = "Pastel1",
                  choices = c("Accent", "Dark 2" = "Dark2", "Paired", 
                              "Pastel 1"= "Pastel1", 
                              "Pastel 2" = "Pastel2", 
                              "Set 1" = "Set1", 
                              "Set 2" = "Set2", "Set 3" = "Set3")),
      actionButton("calculate", label = "Make Figure", 
                   icon("ok", lib = "glyphicon")) #thumbs up icon 
      
      
    ),
    
    #Show the finished plot
    mainPanel(
       tableOutput("taxonomy_table"),
       tableOutput("samples_table"),
       tableOutput("design_table"),
       plotOutput("submit", height = 1000)
       #downloadButton(outputId = "download_data", label = "Download Data")
    )
  )
)

server <- function(input, output) {

  # Download file

  # output$download_data <- downloadHandler(
  # 
  #   filename = function() {
  #     paste0("relative_abundances.", input$filetype)
  #   },
  #   content = function(file) { 
  #     if(input$filetype == "csv"){ 
  #       inFile <- input$taxonomy_file
  #       table <- read.csv(inFile$datapath, header = TRUE)
  #       write_csv(table) 
  #     }
  #     if(input$filetype == "tsv"){ 
  #       inFile <- input$taxonomy_file
  #       table <- read.tsv(inFile$datapath, header = TRUE)
  #       write_tsv(table) 
  #     }
  #   }
  # )

  #renders taxonomy table
  output$taxonomy_table <- renderTable({

    req(input$taxonomy_file)
    req(input$show_taxonomy_table)
    inFile <- input$taxonomy_file
    read.csv(inFile$datapath, header = TRUE)
  })
  
  #renders samples table
  output$samples_table <- renderTable({
    
    req(input$samples_file)
    req(input$show_samples_table)
    inFile.samples <- input$samples_file
    read.table(inFile.samples$datapath, header = TRUE)
  })
  
  #renders design table
  output$design_table <- renderTable({
    
    req(input$design_file)
    req(input$show_design_table)
    inFile.design <- input$design_file
    read.table(inFile.design$datapath, header = TRUE)
  })
  
  #prints structure of taxonomy file
  output$taxonomy <- renderPrint({
    req(input$taxonomy_file)
    str(input$taxonomy_file)
  })
  
  #prints structure of samples file
  output$samples <- renderPrint({
    req(input$samples_file)
    str(input$samples_file)
  })

  #prints structure of design file
  output$design <- renderPrint({
    req(input$design_file)
    str(input$design_file)
  })
  
  #action submit button
  output$submit <- renderPlot({

    req(input$taxonomy_file)
    req(input$calculate) #requires submission
    output$status <- renderPrint({
      print("Submission Received...  Processing Data...")
    })
    
    TIMEPOINT <- input$timepoint
    TITLE <- input$title
    inFile.taxonomy <- input$taxonomy_file
    inFile.samples <- input$samples_file
    inFile.design <- input$design_file
#R CODE#################################################################################
    tax <- read.table(file=inFile.taxonomy$datapath,
                      row.names = 1,
                      header=TRUE,
                      check.names=FALSE,
                      comment.char="") #will become otu.class
    otu.info <- read.table(file=inFile.samples$datapath, header=TRUE, row.names = 2) #will become rm_g
    meta <- read.table(file=inFile.design$datapath, row.names = 1, header =TRUE)
    output$status <- renderPrint({
      print("Files accepted.")
    })
    #make tax and otu.class same dimensions.  cut out all redundant otus
    
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
    otu.class <- data.frame(matrix(NA, nrow = TAXONOMY.LENGTH, ncol = 5))
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
    output$status <- renderPrint({
      print("Taxonomy Information Obtained and Edited.")
    })
    ###CREATE RM_G FILE###
    ###PART 1. CREATE META COLOR FILE###
    taxonomy.color <- data.frame(matrix(NA, nrow = TAXONOMY.LENGTH, ncol = 4)) #4 for OTU, Phylum, Genus, Color
    colnames(taxonomy.color) <- c("OTU", "Phylum", "Genus", "Color")
    taxonomy.color$OTU <- rownames(tax) #add otu numbers
    taxonomy.color$Phylum <- otu.class$Phylum
    taxonomy.color$Genus <- otu.class$Genus
    colors <- as.vector(colorRampPalette(brewer.pal(9, input$color_list))(9))
    colors[10] <- '#ef8f8d' #add a tenth color becaue pastel1 only offers 9
    # colors <- rev(colors)
    output$status <- renderPrint({
      print("Color scheme accepted")
    })

    #define length of vector for the "for loop"
    length.color <- as.numeric(nrow(taxonomy.color))
    
    #ran table(taxonomy.color$Phylum) to get all different phylum present
    #and in what frequency
    #created color palette with as many colors as different phylum present
    #assign colors based on phylum, non-conforming samples are assigned black
    for(i in 1:length.color){
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
    output$status <- renderPrint({
      print("Colors assigned to OTUs")
    })

    ###PART 2. CREATE RELATIVE ABUNDANCE FILE###
    otu.matrix <- as.matrix(otu.info)
    otu.rel <- otu.matrix/rowSums(otu.matrix)
    #otu.rel <- subset(otu.rel, select = -c(label, numOtus))
    otu.rel.t <- t(otu.rel) #transpose
    otu.rel.t <- as.data.frame(otu.rel.t) #make it so we can add back in OTU names without changing to list
    otu.rel.t$OTU <- rownames(otu.rel.t) #add OTUs so we can merge
    rm_g <- merge(taxonomy.color, otu.rel.t, by.x = "OTU", by.y = "OTU")
    ###PART 2. COMPLETE###
    ### rm_g FILE COMPLETED ###
    output$status <- renderPrint({
      print("Merging data tables successful!")
    })

    
    otubar <- as.matrix(subset(rm_g, select =-c(Genus, Color)))
    rownames(otubar) <- otubar[,"Phylum"]
    otubar <- subset(otubar, select = -c(Phylum, OTU))
    bar <- as.data.frame(t(otubar))
    bar$SampleID <- meta$SampleID
    bar$Group <- meta$Group
    
    # barg$topotu <- names(barg)[apply(barg, 1, which.max)] #find the highest otu and add it to topotu
    # barg <- merge(barg, otu.class, by.x=c("topotu"), by.y=c("OTU"))
    col.gen <- as.character(rm_g$Color)
    # bar <- merge(barg, meta, by.x = c("SampleID"), by.y = c("SampleID")) #note: you MUST put barg first! otherwise, it will merge incorrectly
    
    bar_ordered<- bar[order(bar$Group, bar$SampleID),]
    #splits mets and controls
    all<-split(bar_ordered, bar_ordered$Group) 
    met<-all$'MetPN' 
    ctrl<-all$'CtrlPN'
    MET.LENGTH <- as.numeric(nrow(met))
    CTRL.LENGTH <- as.numeric(nrow(ctrl))
    output$status <- renderPrint({
      print("Met and Ctrl tables created")
    })

    # write.table(sums.total, "test.txt", quote=FALSE , sep="\t", col.names=NA)
    ###MAKE MET FILE###
    met <- subset(met, select = -c(SampleID, Group))
    met.t <- t(met)
    met.t <- as.data.frame(met.t)
    met.t$Phylum <- rm_g$Phylum
    met.t <- as.matrix(met.t)
    rownames(met.t) <- met.t[,"Phylum"]
    barmet <- met.t[order(met.t[,"Phylum"]),] #sort matrix by phylum
    barmet <- subset(barmet, select = -c(Phylum)) #might need this...
    class(barmet) <- "numeric" #change matrix to numeric form rather than character
    # barmet: met.t without phylum column
    # met.t: all data
    
    ###phyla only distribution###
    barmet <- as.data.frame(barmet)
    barmet <- t(barmet)
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
    output$status <- renderPrint({
      print("Met sums calculated")
    })
    ###MAKE CTRL FILE###
    ctrl <- subset(ctrl, select = -c(SampleID, Group))
    ctrl.t <- t(ctrl)
    ctrl.t <- as.data.frame(ctrl.t)
    ctrl.t$Phylum <- rm_g$Phylum
    ctrl.t <- as.matrix(ctrl.t)
    rownames(ctrl.t) <- ctrl.t[,"Phylum"]
    barctrl <- ctrl.t[order(ctrl.t[,"Phylum"]),] #sort matrix by phylum
    barctrl <- subset(barctrl, select = -c(Phylum)) #might need this...
    class(barctrl) <- "numeric" #change matrix to numeric form rather than character
    # barctrl: ctrl.t without phylum column
    # ctrl.t: all data
    
    ###phyla only distribution###
    barctrl <- as.data.frame(barctrl)
    barctrl <- t(barctrl)
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
    output$status <- renderPrint({
      print("Ctrl sums calculated")
    })
    sums.total <- cbind(sums.t, sums.ctrl.t)
    output$status <- renderPrint({
      print("Sums Tables Combined")
    })
    colnames(sums.total) <- c(rep("Met PN", MET.LENGTH), rep("Ctrl PN", CTRL.LENGTH))
    # graphing both sets:
    par(mfrow=c(1,1)) #creates a plot of two columns and one row
    par(mar=c(3.3,3,2,1))
    par(xpd=T)
    output$status <- renderPrint({
      print("Producing Graph!")
    })
    plots <- barplot(sums.total, 
                     las=2, 
                     main=TITLE, 
                     ylab="Relative Abundance", 
                     cex.names=.8, 
                     ylim=c(0,1), 
                     col=colors, 
                     xlim=c(0,40),
                     space=c(rep(0.2, MET.LENGTH), 1.5, rep(0.2, CTRL.LENGTH-1)))
    legend.x <- MET.LENGTH + CTRL.LENGTH + 5
    legend <- legend(legend.x, 1,
                     legend=rownames(sums.t),
                     col=colors,
                     fill=colors,
                     cex=1, 
                     bty="n", 
                     ncol=1)
    output$status <- renderPrint({
      print("Done.")
    })
    return(plots)
##############################################################################################
  })
}

shinyApp(ui = ui, server = server)

