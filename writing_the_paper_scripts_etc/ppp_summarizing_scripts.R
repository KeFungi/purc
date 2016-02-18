#### Script to summarize the pacbio data for the PPP manuscript

### Folder organization scheme:
# For these functions to work, the data have to be organized in a particular way 
# (organized [poorly] with the "mainDirectory" setting below). Within that mainDirectory the analyses are organized
# by run, e.g., /R2a, /R2b, /R3a, etc.

# The rounds refer to the three different pacbio submissions that we're analysing (R2, R3, and R4)
# The regimes refer to the clustering settings that we're comparing (calling them a, b, c, etc):

library(ape)
#library(logspline)
#library(gdata)
#library(xlsx)

figureDirectory <- "/Users/carlrothfels/Box Sync/Cystopteridaceae_projects/pacbio_lowcopy_cystopteridaceae/_aa_ppp_manuscript/figures/"
mainDirectory <- "/Users/carlrothfels/Desktop/PURC_analysesForPaper/"


#### Producing histograms of the number of expected errors per sequence
# Using the data from each run that are >600bp
# Example usearch command to get the ".._uncleaned.fa" files:
# usearch7.0.1090 -fastq_filter R2_reads_of_insert.fastq -fastaout R2_uncleaned.fa -eeout -fastq_minlen 600

producePlot_ees <- function(){
  setwd(paste(as.character(mainDirectory), "rawdata/R2", sep = "")) # Moving into a subdirectory of the main directory
  R2 <- read.dna("R2_uncleaned.fa", format = "fasta")
  R2_labels <- names(R2)
  length(R2_labels)
  R2_ees <- as.numeric(sub(".*;ee=(.*);", "\\1", R2_labels))
  
  setwd(paste(as.character(mainDirectory), "rawdata/R3", sep = ""))
  R3 <- read.dna("R3_uncleaned.fa", format = "fasta")
  R3_labels <- names(R3)
  length(R3_labels)
  R3_ees <- as.numeric(sub(".*;ee=(.*);", "\\1", R3_labels))
  
  setwd(paste(as.character(mainDirectory), "rawdata/R4", sep = ""))
  R4 <- read.dna("R4_uncleaned.fa", format = "fasta")
  R4_labels <- names(R4)
  length(R4_labels)
  R4_ees <- as.numeric(sub(".*;ee=(.*);", "\\1", R4_labels))
  
  # length(R4_ees[R4_ees < 8]) # The numbers generated here are close to, but 
  # not the same as, those generated from the usearch quality control runs, somehow
  
  setwd(figureDirectory) # Will save the figure to this directory
  pdf("ee_plots_test1.pdf", h = 2, w = 8)
  breaks = seq(0, 600, by=0.2)
  xlim=c(0,10)
  ylim=c(0,6000)
  par(mfcol = c(1,3),  mar = c(2,2,2,0.5))
  
  # The [R2_ees <10] bit is because it seems that uclust starts rounding the
  # ee numbers starting at 10, so there is a small artifactual peak at 10.0
  h2 <- hist(R2_ees[R2_ees < 10], breaks = breaks, plot=FALSE)
  h3 <- hist(R3_ees[R3_ees < 10], breaks = breaks, plot=FALSE)
  h4 <- hist(R4_ees[R4_ees < 10], breaks = breaks, plot=FALSE)
  
  plot(h2, xlim=xlim, ylim=ylim, main="")
  plot(h3, xlim=xlim, ylim=ylim, main="")
  plot(h4, xlim=xlim, ylim=ylim, main="")
  
  dev.off()
} # End of function


#### Producing the table of coverage/sequence (i.e., per cluster)

## Function to pull the allele count and coverages (size=) from an alignment
summarize_coverages <- function(directory, infile){
  setwd(paste(as.character(mainDirectory), "6_runsWchimeraCountsOutput/", as.character(directory), sep = ""))
  file <- read.dna(infile, format = "fasta")
  seqnames <- names(file)
  coverages <- as.numeric(sub(".*size=(.*)", "\\1", seqnames))
  clusters <- length(coverages)
  coverageMean <- mean(coverages)
  coverageSD <- sd(coverages)
  output <- list(clusters, coverageMean, coverageSD, coverages)
  return(output)
} # End of function

## Function that calls summarize_convergences for all the analyses and puts the results together in a table
produceTable_coverages <- function(regimes, loci){
  coverage_table <- matrix(nrow=length(regimes)*length(rounds), ncol=length(loci))
  colnames(coverage_table) <- loci
  
  # An unnecessarily (??) complicated way to get rownames in the form of "regime_round"
  rnames <- list() 
  count = 0
  for (re in regimes){
    for (ro in rounds){
      count  <- count +1
      rnames[[count]] <- paste(re, ro, sep="_")
    }
  }
  rownames(coverage_table) <- rnames
  
  rcount = 0
  for (regime in regimes){
    for (round in rounds){
      ccount = 0
      rcount = rcount + 1 # What row are we in
      #directory <- paste(round, "/", round, regime, sep="") # Getting directories in the form of, e.g., "R2/R2A"
      directory <- paste(round, regime, sep="") # Getting directories in the form of, e.g., "R2A"
      for (locus in loci){
        ccount = ccount +1
        values <- summarize_coverages(directory, paste(locus, "clustered_renamed.fa", sep="_"))
        cell <- paste(values[[1]], " (", round(values[[2]], digits=1), " +/- ", round(values[[3]], digits=1), ")", sep="")
        coverage_table[rcount,ccount] <- cell
      }
    } 
  }
  return(coverage_table)
} # End of function


#### Producing plots (histograms) of coverage values

## Function to create (but not plot) a histogram of coverage depths/allele
# May want to do this for the "raw" and "corrected" set of alleles? To show that the bad ones have lower coverage?
# Not currently using this -- the curves ones, below, seem to be better
producePlot_coverages <- function(directory, infile){
  values <- summarize_coverages(directory, infile)
  coverages <- values[[4]] # summarize_coverages() returns four things -- clusters, coverageMean, coverageSD, coverages -- and we only want the last one here
  breaks  <- seq(0, 200, by=5) #Starting the breaks at 4 because that's the minumum size allowed under my PURC settings for these runs
  # Need to tweak the xlim etc based on the result. Currently all the coverages >200 are ignored
  h <- hist(coverages[coverages<201], breaks =  breaks, plot = FALSE)
  toReturn <- list(h, coverages) # Need to return coverages here inorder to be able to call abline in the other function
  return(toReturn) 
} # End of function

## Function like above, but the breaks in this one work better for small numbers of plots
producePlot_coverages_old <- function(directory, infile){
  values <- summarize_coverages(directory, infile)
  coverages <- values[[4]]
  breaks  <- seq(4, max(coverages)+5, by=2) #Starting the breaks at 4 because that's the minumum size allowed under my PURC settings for these runs
  h <- hist(coverages, breaks =  breaks, plot = FALSE)
  return(h)
} # End of function


## Function producing initial curves (and/or histograms) matching the counts rather than the densities
# These two fxns are working great, except that the curves go down toward the origin, even though there are no data with coverage >4. Haven't been able to figure out how to make it a truncated density.
producePlot_coveragesWcurve <- function(directory, infile, colour = "blue", add, title, binWidth = binWidth){
  values <- summarize_coverages(directory, infile)
  coverages <- values[[4]] 
  newDensity <- density(coverages) #adjust = 1.2
  newDensity$y <- newDensity$y * length(coverages) * binWidth # The Yvalues of the density produce a curve that sums to one
  # So multipling them by the number of alleles and the histogram binwidth gets a curve that matches the coverage counts
  
  hist(coverages[coverages<121], breaks = seq(4,120, by=binWidth), xlim = c(0, 120), ylim = c(0,35), ylab= "Number of Alleles", xlab= "Coverage", main = title, col = "lightgrey")
  lines(newDensity, col=colour, lwd=2, xlim= c(4,120))
  
  # Either do the hist + lines, above, or just the plot, below (for a figure without the underlying histogram)
  # plot(newDensity, col=colour, lwd=2, xlim = c(4, 120), ylim = c(0,80), ylab= "Number of Alleles", xlab= "Coverage", main = title) 
} # End of function

## Function that adds density curves to a preexisting plot (i.e., from producePlot_coveragesWcurve)
addLines_coveragesWcurve <- function(directory, infile, colour = "blue", binWidth = binWidth){
  values <- summarize_coverages(directory, infile)
  coverages <- values[[4]]
  newDensity <- density(coverages) #adjust = 1.2
  newDensity$y <- newDensity$y * length(coverages) * binWidth # The Yvalues of the density produce a curve that sums to one
  lines(newDensity, col=colour, lwd=2, xlim = c(4,120)) 
} # End of function


## Master function that goes through each of the regimes, loci, and rounds and plots density curves of the coverages (reads per "allele")
#Those from the regimes for a given locus are plotted on top of each other, in different colours
produceFigure_coveragesWcurve <- function(){ #regimes, loci, rounds
  regimes <- list("a", "c", "e") #"b", "d",
  loci <- list("APP", "IBR", "PGI", "GAP") #"GAP"
  rounds <- list("R2") #, "R3", "R4" not done yet
  colours  <-  c("red", "orange", "blue", "purple", "pink") # A master list of colours
  colours  <- colours[1:length(regimes)] # Getting one colour for each of the regimes
  binWidth = 3
  
  # mfrow=c(nrows, ncols) # fills in by row
  #dev.off() # Resetting the plot device
  makepdf = TRUE
  if (makepdf){
    setwd(figureDirectory)
    pdf("coverage_plots_curves2.pdf", h = 4, w = 7) #, h = 4, w = 8
    setwd(mainDirectory)
  }
  par(mfrow=c(length(rounds), length(loci)))
  
  for (locus in loci){
    for (round in rounds){
      regimeNum  <- 0
      for (regime in regimes){
        regimeNum  <- regimeNum + 1
        directory <- paste(round, regime, sep="")
        infile <- paste(locus, "clustered_renamed.fa", sep="_")
        if (regimeNum == 1){
          producePlot_coveragesWcurve(directory, infile, colour = colours[regimeNum], add=FALSE, title=locus, binWidth = binWidth)
          #abline(v=4) # To show the cut-off (only kept alleles with coverage >4)
        } else {
          addLines_coveragesWcurve(directory, infile, colour = colours[regimeNum], binWidth = binWidth)
        }
      }
    }
  }
  dev.off()
} # End of function


## workspace functions -- fooling around
fooling_tryingToGetCurvesOnCountData <- function(){
  x <- 1:10
  y <- c(2,4,6,8,7,12,14,16,18,20)
  lo <- loess(y~x)
  plot(x,y)
  lines(predict(lo), col='red', lwd=2)
  
  dev.off()
  h = hist(y, breaks = seq(0,20, by=2))
  tocurve <- loess(h$counts ~ h$mids)
  #plot(h$mids, h$counts)
  #lines(predict(tocurve), col="red") # This gets thrown off by the by=2 in breaks and produces a line that's horizontally compressed
  
  lines(seq(1,20,by=2), tocurve$fitted) #This works, awkwardly
  
  # Or this -- has _some_ potential
  round <- "R2"
  regime <- "a"
  locus <- "APP"
  directory <- paste(round, regime, sep="")
  infile <- paste(locus, "clustered_renamed.fa", sep="_")
  values <- summarize_coverages(directory, infile)[[4]]
  h <- hist(values, breaks = seq(0,81, by=3))
  
  smoo <- spline(h$breaks[1:length(h$breaks)-1],h$counts)
  lines(smoo$x[-1],smoo$y[-1],col='green') # The -1s are an attempt to get rid of the first elements, and thus not have the curve go to the origin
}

#Function to test the producePlot_coverages() function
exampleCovPlot <- function(){
  round <- "R2"
  regime <- "a"
  locus <- "APP"
  directory <- paste(round, regime, sep="")
  infile <- paste(locus, "clustered_renamed.fa", sep="_")
  plot(producePlot_coverages(directory, infile)[[1]], main ="", xlab ="", ylab ="", col = "lightblue") 
} # End of function
exampleCovPlot()


## Function that calls producePlot_coverages for all the analyses and puts the
# results together in a figure -- analgous to the produceTable_coverages function
produceFigure_coverages <- function(regimes, loci, rounds){
  # mfrow=c(nrows, ncols) # fills in by row
  par(mfrow=c(length(loci), length(rounds)*length(regimes)))
  
  for (locus in loci){
    for (round in rounds){
      for (regime in regimes){
        #directory <- paste(round, "/", round, regime, sep="")
        directory <- paste(round, regime, sep="")
        infile <- paste(locus, "clustered_renamed.fa", sep="_")
        results <- producePlot_coverages(directory, infile)
        toPlot <- results[[1]] #PLoting the histogram, rather than the second element which is the coverages counts
        plot(toPlot, main = paste(locus, round, regime, sep="_"), xlab="", ylab = "", ylim=c(0,20), xlim=c(0,200)) 
        abline(v=median(results[[2]])) # The second element of toPlot has the lost if coverage values
      }
    }
  }
} # End of function


### Functions for dealing with the chimeras

makeChimeraPlot <- function(){
  chimeraDF <- read.csv("/Users/carlrothfels/Box Sync/R_Python_etal/git_repositories/ppp_repo/writing_the_paper/chimeraCounts.csv", header=FALSE, stringsAsFactors = FALSE)
  chimeraDF$V7 = NULL # An extra column gets in somehow -- erasing it here
  colnames(chimeraDF) <- c("regime", "locus","CK1", "CK2", "CK3", "CK4") #"CK" for "chimera kill"
  regimes <- sort(unique(chimeraDF$regime))
  loci <- sort(unique(chimeraDF$locus))
  
  par(mfrow=c(length(loci), 1), mar=c(1,2,1,2), oma = c(5,5,0,0)) #mar - margins of indivual plots, b,L,t,R
  for (thislocus in loci){
    data <- subset(chimeraDF, locus == thislocus)
    data <- data[-c(1,2)] # Getting rid of the first two columns
    data <- data.matrix(data) # Changing it from a dataframe to a matrix
    barplot(data, beside=TRUE, ylim= c(0,200), main=thislocus)
  }
  title(xlab = "Clustering Regime for Each Round of Chimera Detection",
        ylab = "Number of Chimeras Detected", outer = TRUE, line = 3)
} # End of function

oldmakeChimeraPlot <- function(){
  colours <- c("red", "blue", "orange", "red", "blue", "orange")
  for (thislocus in loci){
    data <- subset(chimeraDF, locus == thislocus)
    regimeCount  <- 0
    for (Regime in regimes){
      regimeCount = regimeCount + 1
      row  <- subset(data, regime == Regime)
      cat("Regime is ", Regime)
      if (Regime == "R2a"){
        #smoothed <- loess(as.integer(row[3:6]) ~ seq(1,4))
        plot(seq(1,4), row[3:6], ylim = c(0,150), type = "l", xlab ="") #, add = TRUE
        #lines(predict(smoothed))
        lines(lowess(seq(1,4), row[3:6]), col = colours[regimeCount], lwd =2)
      } else{
        if (Regime == "R2c" | Regime == "R2e"){
          lines(lowess(seq(1,4), row[3:6]), col = colours[regimeCount], lwd =2)
        }else{ # Want these ones to have dotted lines
          lines(lowess(seq(1,4), row[3:6]), col = colours[regimeCount], lty = 2, lwd = 2)
        }
        
      }
    }
  }
} #End of function

#### Calls

# To plot the ee values:
producePlot_ees()

# To summarize the coverage values for all the runs
regimes <- list("a", "c", "e") #"b", "d",
stringentRegimes <- paste("Str_", regimes, sep="") # making list of e.g., "Str_a" regimes
#regimes = c(regimes, stringentRegimes)
loci <- list("APP", "GAP", "IBR", "PGI")
rounds <- list("R2") #, "R3", "R4" not done yet

table <- produceTable_coverages(regimes, loci)

# To produce a plot of coverage histograms
setwd(figureDirectory) # Will save the figure to this directory
pdf("coverage_plots_test_d.pdf") #, h = 7, w = 8
produceFigure_coverages(regimes, loci, rounds)
dev.off()

# To produce the coverage plots, one for each locus, with coloured lines showing the performance of the different regimes
produceFigure_coveragesWcurve()

# To produce a barplot of the chimera counts
makeChimeraPlot()
