#### Script to summarize the pacbio data for the PPP manuscript

#Not sure why Rstudio won't let me open this from within the project
#I have to right-click it and "open with Rstudio"
# (maybe this is fixed now)

### Folder organization scheme:
# For these functions to work, the data have to be organized in a particular way 
# (organized [poorly] with the "mainDirectory" setting below)
## The uncleaned data (for plotting numbers of expected errors) are at:
# mainDirectory/R2/1_calcultngEEs_noQC
# mainDirectory/R3/1_calcultngEEs_noQC   etc
## Then the PURC outputs are in mainDirectory/round/regime
# e.g., mainDirectory/R2/R2a
#       mainDirectory/R2/R2b   etc

# The rounds refer to the three different pacbio submissions that we're analysing (R2, R3, and R4)
# The regimes refer to the clustering settings that we're comparing (calling them A, B, C):
# A: 0.999, 0.997, 0.990. B: 0.997, 0.995, 0.990. C: 0.990, 0.980, 0.970
# (for all regimes, all clusters are retained after clusterings 1 and 2, and only 
# clusters >3sequences in size after clustering 3.)


library(ape)

figureDirectory <- "/Users/carlrothfels/Box Sync/Cystopteridaceae_projects/pacbio_lowcopy_cystopteridaceae/_aa_ppp_manuscript/figures/"
mainDirectory <- "/Users/carlrothfels/Box Sync/Cystopteridaceae_projects/pacbio_lowcopy_cystopteridaceae/aa_data_and_analysis/201509_analysesformanuscript/"


#### Producing histograms of the number of expected errors per sequence
# Using the data from each run that are >600bp
# Example usearch command to get the ".._uncleaned.fa" files:
# usearch7.0.1090 -fastq_filter R2_reads_of_insert.fastq -fastaout R2_uncleaned.fa -eeout -fastq_minlen 600

producePlot_ees <- function(){
  setwd(paste(as.character(mainDirectory), "R2/1_calcultngEEs_noQC", sep = "")) # Moving into a subdirectory of the main directory
  R2 <- read.dna("R2_uncleaned.fa", format = "fasta")
  R2_labels <- names(R2)
  length(R2_labels)
  R2_ees <- as.numeric(sub(".*;ee=(.*);", "\\1", R2_labels))
  
  setwd(paste(as.character(mainDirectory), "R3/1_calcultngEEs_noQC", sep = ""))
  R3 <- read.dna("R3_uncleaned.fa", format = "fasta")
  R3_labels <- names(R3)
  length(R3_labels)
  R3_ees <- as.numeric(sub(".*;ee=(.*);", "\\1", R3_labels))
  
  setwd(paste(as.character(mainDirectory), "R4/1_calcultngEEs_noQC", sep = ""))
  R4 <- read.dna("R4_uncleaned.fa", format = "fasta")
  R4_labels <- names(R4)
  length(R4_labels)
  R4_ees <- as.numeric(sub(".*;ee=(.*);", "\\1", R4_labels))
  
  # length(R4_ees[R4_ees < 8]) # The numbers generated here are close to, but 
  # not the same as, those generated from the usearch quality control runs, somehow
  
  setwd(figureDirectory) # Will save the figure to this directory
  pdf("ee_plots_test.pdf", h = 2, w = 8)
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

# Function to pull the allele count and coverages (size=) from an alignment
summarize_coverages <- function(directory, infile){
  setwd(paste(as.character(mainDirectory), as.character(directory), sep = ""))
  file <- read.dna(infile, format = "fasta")
  seqnames <- names(file)
  coverages <- as.numeric(sub(".*size=(.*)", "\\1", seqnames))
  clusters <- length(coverages)
  coverageMean <- mean(coverages)
  coverageSD <- sd(coverages)
  output <- list(clusters, coverageMean, coverageSD, coverages)
  return(output)
} # End of function

# function that calls summarize_convergences for all the analyses and puts the
# results together in a table
produceTable_coverages <- function(regimes, loci){
  coverage_table <- matrix(nrow=length(regimes)*length(rounds), ncol=length(loci))
  colnames(coverage_table) <- loci
  
  rnames <- list() # An unnecessarily (??) complicated way to get rownames in the form of "regime_round"
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
      directory <- paste(round, "/", round, regime, sep="") # Getting directories in the form of, e.g., "R2/R2A"
      for (locus in loci){
        ccount = ccount +1
        ## need to change this values so that the directory passed isn't just regime, but instead paste(regime, round)
        values <- summarize_coverages(directory, paste("_", locus, "_clustered_renamed.fa", sep=""))
        cell <- paste(values[[1]], " (", round(values[[2]], digits=1), " +/- ", round(values[[3]], digits=1), ")", sep="")
        coverage_table[rcount,ccount] <- cell
      }
    } 
  }
  return(coverage_table)
} # End of function


#### Producing plots (histograms) of coverage values

# Function to produce a plot of coverage depths/allele
# May want to do this for the "raw" and "corrected" set of alleles? To show that the bad ones have lower coverage?
producePlot_coverages <- function(directory, infile){
  values <- summarize_coverages(directory, infile)
  coverages <- values[[4]]
  breaks  <- seq(0, 200, by=5) #Starting the breaks at 4 because that's the minumum size allowed under my PURC settings for these runs
  # Need to tweak the xlim etc based on the result. Currently all the coverages >200 are ignored
  h <- hist(coverages[coverages<201], breaks =  breaks, plot = FALSE)
  toReturn <- list(h, coverages) # Need to return coverages here inorder to be able to call abline in the other function
  return(toReturn) 
} # End of function

# The breaks in this one work better for small numbers of plots
producePlot_coverages_old <- function(directory, infile){
  values <- summarize_coverages(directory, infile)
  coverages <- values[[4]]
  breaks  <- seq(4, max(coverages)+5, by=2) #Starting the breaks at 4 because that's the minumum size allowed under my PURC settings for these runs
  h <- hist(coverages, breaks =  breaks, plot = FALSE)
  return(h)
} # End of function

#Function to test the producePlot_coverages() function
exampleCovPlot <- function(){
  round <- "R2"
  regime <- "A"
  locus <- "APP"
  directory <- paste(round, "/", round, regime, sep="")
  infile <- paste("_", locus, "_clustered_renamed.fa", sep="")
  plot(producePlot_coverages(directory, infile)[[1]], main ="", xlab ="", ylab ="", col = "lightblue") 
} # End of function
exampleCovPlot()


# function that calls producePlot_coverages for all the analyses and puts the
# results together in a figure -- analgous to the produceTable_coverages function
produceFigure_coverages <- function(regimes, loci, rounds){
  # mfrow=c(nrows, ncols) # fills in by row
  par(mfrow=c(length(loci), length(rounds)*length(regimes)))
  
  for (locus in loci){
    for (round in rounds){
      for (regime in regimes){
        directory <- paste(round, "/", round, regime, sep="")
        infile <- paste("_", locus, "_clustered_renamed.fa", sep="")
        results <- producePlot_coverages(directory, infile)
        toPlot <- results[[1]] #PLoting the histogram, rather than the second element which is the coverages counts
        plot(toPlot, main = paste(locus, round, regime, sep="_"), xlab="", ylab = "", ylim=c(0,20), xlim=c(0,200)) 
        abline(v=median(results[[2]])) # The second element of toPlot has the lost if coverage values
      }
    }
  }
} # End of function



#### Calls

# To plot the ee values:
producePlot_ees()

# To summarize the coverage values for all the runs
regimes <- list("A", "B", "C")
loci <- list("APP", "GAP", "IBR", "PGI")
rounds <- list("R2") #, "R3", "R4" not done yet

table <- produceTable_coverages(regimes, loci)

# To produce a plot of coverage histograms
setwd(figureDirectory) # Will save the figure to this directory
pdf("coverage_plots_test.pdf", h = 7, w = 8)
produceFigure_coverages(regimes, loci, rounds)
dev.off()

