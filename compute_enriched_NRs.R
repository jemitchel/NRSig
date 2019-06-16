library(annotate)
library(jsonlite)
library("hgu133plus2.db")
library(ggplot2)

# gets gene symbol annotations for all probes
x <- hgu133plus2SYMBOL
mapped_probes <- mappedkeys(x)
glist <- as.list(x[mapped_probes]) #list of symbols with only probes that map to genes

# opens the tf-target probe data source 
allTargets <- fromJSON("C:/Users/jonat/Documents/R/NRSig-app/data/NR_target_probes.json", flatten=TRUE) #use this one for testing
# allTargets <- fromJSON("./data/NR_target_probes.json", flatten=TRUE)

# loads dataframe of serum starved MCF7 samples
data <- readRDS(file = "C:/Users/jonat/Documents/R/NRSig-app/data/serum-starved.rds") #use this for testing
# data <- readRDS(file = "./data/serum-starved.rds") 

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
  mappings <- list()
  for (i in 1:nrow(probeZ)) {
    ID <- rownames(probeZ)[i] 
    if (ID %in% mapped_probes) { #if the probe maps to a gene symbol, continue
      sym <- glist[[ID]]
      score <- probeZ[i,1]
      if (sym %in% names(zForGenes)) {
        zForGenes[[sym]] <- c(zForGenes[[sym]],score)
        mappings[[sym]] <- c(mappings[[sym]],ID)
      } else {
        zForGenes[[sym]] <- c(score)
        mappings[[sym]] <- c(ID)
      }
    }
  }
  
  newz <- list()
  newz_mappings <- list()
  for (i in 1:length(zForGenes)) {
    ndxMax <- which.max(abs(zForGenes[[i]]))
    newz[[names(zForGenes)[i]]] <- zForGenes[[i]][ndxMax] #not that this is not abs value
    newz_mappings[[names(mappings)[i]]] <- mappings[[i]][ndxMax]
  }

  return(list(newz,newz_mappings))
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

MakeBoxPlots <- function(testData,gene_z,mappings,NR) {

  # preps serum-starved data for plotting
  distrib <- data[unlist(mappings),]
  genelst <- rep(glist[rownames(distrib)], each = ncol(distrib))
  distrib <- c(t(distrib)) # linearizes data
  distrib <- cbind(distrib,genelst)
  distrib <- data.frame(distrib)
  colnames(distrib) <- c("expr","gene")
  distrib <- as.data.frame(lapply(distrib, unlist))
  
  # preps user input data and results for plotting
  ordered_z <- as.list(gene_z[order(abs(unlist(gene_z)),decreasing = FALSE)]) # orders the list
  inData <- testData[unlist(mappings),,drop=FALSE]
  genelst <- glist[rownames(inData)]
  inData <- cbind(inData,unlist(genelst))
  colnames(inData) <- c("expr","gene")
  inData <- cbind(inData,unlist(gene_z[as.character(inData$gene)]))
  colnames(inData)[3] <- "zsc"
  inData <- inData[match(names(ordered_z), inData$gene),]
  inData$gene <- factor(inData$gene, levels = inData$gene) # orders the factor for plot
  inData$zsc <- round(inData$zsc,digits = 2)
   
  # makes the plot
  p <- ggplot(data = inData, aes(x=as.factor(gene), y=expr)) +
    geom_point(color='red', size=3) +
    geom_text(data = inData, aes(x=as.factor(gene), y=20, label=zsc), vjust=0) +
    geom_boxplot(data = distrib, aes(x=as.factor(gene), y=expr), width=.3) +
    geom_point(data = inData, aes(x=as.factor(gene), y=expr), color='red', size=3) +
    coord_flip() +
    labs(title=paste(NR, " Target Expression", sep=""), 
         subtitle = "Z-Score", x="Target Gene", y="fRMA Expression") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = .975)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  return(p)
}


# This is the main function
CalcEnrich <- function(testSamples){
  
  # keeps track of runtime for progress bar
  runtimes <- c()
  num_completed <- 0
  
  # averages replicates if there are multiple
  test_sample <- data.frame(Mean=rowMeans(testSamples,na.rm=TRUE))
  
  # compute zscores for all probes (test sample versus serum-starved prior distributions)
  zscores <- (test_sample-av)/std
  
  # creates df to store enrichment p-values
  finPvals <- data.frame(matrix(ncol = length(NRs), nrow = 1))
  colnames(finPvals) <- NRs
  
  finNR <- list() # holds z-scores for target genes of each NR
  
  for (d in 1:length(NRs)) { # for each nuclear receptor...
    
    # gets a list of target probes for the specified NR
    if (NRs[d] %in% names(allTargets)) {
      target_probes <- allTargets[[NRs[d]]]
    } else {
      next # if NR had less than 15 target genes go to next NR
    }
    
    start_time <- Sys.time() # starts timer
    
    # select just the probes that are targets of the current NR being evaluated
    test_sample_nr <- zscores[target_probes,,drop=F]
    
    # select all other probes that are not listed as targets of the NR
    test_sample_rest <- zscores[!row.names(zscores)%in%target_probes,,drop=F]
    
    # gets differentially expressed genes from probe z scores
    nrGeneLevelZ <- DiffProbes2Genes(test_sample_nr)
    nrGeneLevelZ_mappings <- nrGeneLevelZ[[2]]
    nrGeneLevelZ <- nrGeneLevelZ[[1]]
    restGeneLevelZ <- DiffProbes2Genes(test_sample_rest)
    restGeneLevelZ_mappings <- restGeneLevelZ[[2]]
    restGeneLevelZ <- restGeneLevelZ[[1]]
    
    # generate boxplots of the data
    outPlot <- MakeBoxPlots(test_sample,nrGeneLevelZ,nrGeneLevelZ_mappings,NRs[d])
    
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
    finNR[[NRs[d]]] <- list(outPlot,length(nrGeneLevelZ))
    finPvals[1,NRs[d]] <- result[["p.value"]]

    # calculates remaining compute time
    end_time <- Sys.time()
    runtime <- end_time - start_time
    runtimes <- c(runtimes,runtime)
    num_completed <- num_completed + 1
    time_left <- mean(runtimes) * (15 - num_completed) / 60
    minutes_left <- floor(time_left)
    seconds_left <- round((time_left - minutes_left) * 60)
    
    incProgress(amount = 1/15,
                detail = paste("       ",minutes_left, " minutes ",seconds_left," seconds remaining"))
    
  }
  
  # removes NA values
  finPvals <- finPvals[,!apply(is.na(finPvals), 2, any)]
  
  # appends adjusted p-values
  padj <- p.adjust(finPvals[1,],method = "fdr")
  finPvals <- rbind(finPvals,padj)
  row.names(finPvals) <- c("raw p-value","adj p-value")
  
  # transposes and sorts dataframe
  finPvals <- t(finPvals) 
  finPvals <- finPvals[order(finPvals[,2]),]
  finNR <- as.list(finNR[order(finPvals[,2])]) # orders the list of plots
  
  # adds gene symbol annotations to differentially expressed probes list
  syms <- as.character(glist[rownames(zscores)])
  syms <- replace(syms, syms=="NULL", "")
  zscores <- cbind(zscores,syms)
  # sorts differentially expressed probes list
  zscores <- zscores[order(abs(unlist(zscores[,1])), decreasing=TRUE),]
  
  # makes probe ID to first column and names columns (for both downloadable tablesS)
  zscores <- cbind(Row.Names = rownames(zscores), zscores)
  colnames(zscores) <- c("Probeset ID","Z-Score","Gene Symbol")
  testSamples <- cbind(Row.Names = rownames(testSamples), testSamples)
  colnames(testSamples)[1] <- c("Probeset ID")
  

  return(list(finNR,finPvals,zscores,testSamples))
}







