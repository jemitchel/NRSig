library(annotate)
library(jsonlite)
library("hgu133plus2.db")

# gets gene symbol annotations for all probes
x <- hgu133plus2SYMBOL
mapped_probes <- mappedkeys(x)
glist <- as.list(x[mapped_probes]) #list of symbols with only probes that map to genes

# opens the tf-target probe data source 
allTargets <- fromJSON("C:/Users/jonat/Documents/R/NRSig-app/data/NR_target_probes.json", flatten=TRUE) #use this one for testing
# allTargets <- fromJSON("NR_target_probes.json", flatten=TRUE) #use this one for real

# loads dataframe of serum starved MCF7 samples
data <- readRDS(file = "C:/Users/jonat/Documents/R/NRSig-app/data/serum-starved.rds") #use this for testing
# data <- readRDS(file = "data/serum-starved.rds") #use this one for real

# computes mean and standard deviation for serum starved data
data_av_std <- data.frame(Mean=rowMeans(data,na.rm=TRUE),Std=apply(data,1,sd,na.rm=TRUE))
av <- data_av_std[,1]
std <- data_av_std[,2]

# all 48 NRs (some won't be tested since have < 15 target genes listed in TTRUST)
NRs <- c("ESR1", "AR", "NR3C1", "NR3C2", "PGR",
         "NR4A1", "NR4A2", "NR4A3", "NR5A1", "NR5A2", "NR6A1",
         "NR0B1", "NR0B2", "THRA", "THRB", "RARA", "RARB", "RARG",
         "PPARA", "PPARD", "PPARG", "NR1D1", "NR1D2", "RORA", "RORB",
         "RORC", "NR1H4", "NR1H3", "NR1H2", "VDR", "NR1I2",
         "NR1I3", "HNF4A", "HNF4G", "RXRA", "RXRB", "RXRG", "NR2C1",
         "NR2C2", "NR2E1", "NR2E3", "NR2F1", "NR2F2", "NR2F6",
         "ESR2", "ESRRA", "ESRRB", "ESRRG")

# This function converts z scores from probe level -> gene level
DiffProbes2Genes <- function(probeZ) {
  zForGenes <- list()
  for (i in 1:nrow(probeZ)) {
    ID <- rownames(probeZ)[i] 
    if (ID %in% mapped_probes) { #if the probe maps to a gene symbol, continue
      sym <- glist[[ID]]
      score <- probeZ[i,1]
      if (sym %in% names(zForGenes)) {
        zForGenes[[sym]] <- c(zForGenes[[sym]],score)
      } else {
        zForGenes[[sym]] <- c(score)
      }
    }
  }
  
  newz <- list()
  for (i in 1:length(zForGenes)) {
    ndxMax <- which.max(abs(zForGenes[[i]]))
    newz[[names(zForGenes)[i]]] <- zForGenes[[i]][ndxMax] #not that this is not abs value
  }

  return(newz)
}

# This function counts number of genes that are above the z-score threshold
GetNumAboveThreshold <- function(zscores) {
  count <- 0
  thresh <- 2 # this is the zscore threshold
  for (i in 1:length(zscores)) {
    if (abs(zscores[[i]]) > thresh) {
      count <- count + 1
    }
  }
  return(count)
}

# This is the main function
CalcEnrich <- function(testSamples){
  
  # creates empty dataframe to store enrichment results (just p-values)
  findb <- data.frame(matrix(ncol = length(NRs), nrow = ncol(testSamples)))
  row.names(findb) <- colnames(testSamples)
  colnames(findb) <- NRs
  
  # iterates through each test sample and nuclear receptor
  finSmpl <- list() #will hold z-scores for target genes of each NR for each sample
  for (p in 1:ncol(testSamples)) { #for each sample...
    finNR <- list() #will hold z-scores for target genes of each NR
    test_sample <- testSamples[,p,drop=F] #get epression values for all probes for the current sample
    
    # compute zscores for all probes (test sample versus serum-starved prior)
    zscores <- (test_sample-av)/std #WILL WANT TO RETURN THESE FOR USER
    
    for (d in 1:length(NRs)) { #for each nuclear receptor...
      # gets a list of target probes for the specified NR
      if (NRs[d] %in% names(allTargets)) {
        target_probes <- allTargets[[NRs[d]]]
      } else {
        next # if NR had less than 15 target genes go to next NR
      }
      
      # select just the probes that are targets of the current NR being evaluated
      test_sample_nr <- zscores[target_probes,,drop=F]
      
      # select all other probes that are not listed as targets of the NR
      test_sample_rest <- zscores[!row.names(zscores)%in%target_probes,,drop=F]
      
      # gets differentially expressed genes from probe z scores
      nrGeneLevelZ <- DiffProbes2Genes(test_sample_nr)
      restGeneLevelZ <- DiffProbes2Genes(test_sample_rest)
      
      # figure out how many target and non-targets are dysregulated
      nrAboveThresh <- GetNumAboveThreshold(nrGeneLevelZ)
      restAboveThresh <- GetNumAboveThreshold(restGeneLevelZ)
      
      # design contingency table
      aboveThresh <- c(nrAboveThresh,restAboveThresh)
      belowThresh <- c(length(nrGeneLevelZ)-nrAboveThresh,length(restGeneLevelZ)-restAboveThresh)
      cTable <- cbind(aboveThresh,belowThresh)
      
      # computes results
      result <- fisher.test(cTable,alternative = "greater")
      print(NRs[d])
      print(result[["p.value"]])
      finNR[[NRs[d]]] <- nrGeneLevelZ
      findb[colnames(testSamples)[p],NRs[d]] <- result[["p.value"]]
    }
    finSmpl[[colnames(testSamples)[p]]] <- finNR
  }
}







