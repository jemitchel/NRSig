library(annotate)
library(jsonlite)
library("hgu133plus2.db")
library(ggplot2)
library(genefilter)

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


# Computes T-test across rows
getDiffex <- function(ss,test) {
  combined <- as.matrix(data.frame(ss,test))
  
  # # if assume equal variance TRY THIS ON WHOLE DATASET TO SEE WHICH IS BETTER
  # groups <- as.factor(c(rep(0, ncol(ss)), rep(1, ncol(test))))
  # pvals <- rowttests(combined, groups)
  
  # if don't assume equal variances
  print(64+ncol(test))
  tStats <- fastT(combined,1:64,65:(64+ncol(test)),var.equal=FALSE)
  pvals <- 2*pt(-abs(tStats$z),df=ncol(combined)-2)

  # does fdr adjust on pvalues
  padj <- p.adjust(pvals,method = "fdr")
  
  scores <- data.frame('pval' = padj)
  rownames(scores) <- rownames(ss)
  return(scores)
}

# This function converts z scores from probe level -> gene level
DiffProbes2Genes <- function(probeSc,useZsc) {
  scForGenes <- list()
  mappings <- list()
  for (i in 1:nrow(probeSc)) {
    ID <- rownames(probeSc)[i] 
    if (ID %in% mapped_probes) { #if the probe maps to a gene symbol, continue
      sym <- glist[[ID]]
      score <- probeSc[i,1]
      if (sym %in% names(scForGenes)) {
        scForGenes[[sym]] <- c(scForGenes[[sym]],score)
        mappings[[sym]] <- c(mappings[[sym]],ID)
      } else {
        scForGenes[[sym]] <- c(score)
        mappings[[sym]] <- c(ID)
      }
    }
  }
  
  newz <- list()
  newz_mappings <- list()
  if (useZsc) {
    for (i in 1:length(scForGenes)) {
      ndxMax <- which.max(abs(scForGenes[[i]]))
      newz[[names(scForGenes)[i]]] <- scForGenes[[i]][ndxMax] #not that this is not abs value
      newz_mappings[[names(mappings)[i]]] <- mappings[[i]][ndxMax]
    }
  } else {
    for (i in 1:length(scForGenes)) {
      ndxMin <- which.min(scForGenes[[i]])
      newz[[names(scForGenes)[i]]] <- scForGenes[[i]][ndxMin] #not that this is not abs value
      newz_mappings[[names(mappings)[i]]] <- mappings[[i]][ndxMin]
    }
  }
  
  
  return(list(newz,newz_mappings))
}

# This function counts number of genes that are above the z-score threshold
GetNumAboveThreshold <- function(scores,useZsc) {
  count <- 0
  if (useZsc) {
    thresh <- 2 
    for (i in 1:length(scores)) {
      if (abs(scores[[i]]) > thresh) {
        count <- count + 1
      }
    }
  } else{
    thresh <- .05 
    for (i in 1:length(scores)) {
      if (abs(scores[[i]]) < thresh) {
        count <- count + 1
      }
    }
  }
  
  return(count)
}

MakeBoxPlots <- function(testData,gene_z,mappings,NR,useZsc) {
  
  # preps serum-starved data for plotting
  distrib <- data[unlist(mappings),]
  genelst <- rep(glist[rownames(distrib)], each = ncol(distrib))
  distrib <- c(t(distrib)) # linearizes data
  distrib <- cbind(distrib,genelst)
  distrib <- data.frame(distrib)
  colnames(distrib) <- c("expr","gene")
  distrib <- as.data.frame(lapply(distrib, unlist))
  
  # preps user input data and results for plotting
  testData <- data.frame(Mean=rowMeans(testData,na.rm=TRUE)) #gets means
  if (useZsc) {
    ordered_z <- as.list(gene_z[order(abs(unlist(gene_z)),decreasing = FALSE)]) # orders the list
  } else {
    ordered_z <- as.list(gene_z[order(abs(unlist(gene_z)),decreasing = TRUE)]) # orders the list
    
  }
  inData <- testData[unlist(mappings),,drop=FALSE]
  genelst <- glist[rownames(inData)]
  inData <- cbind(inData,unlist(genelst))
  colnames(inData) <- c("expr","gene") # this is where it breaks because now have more columns
  inData <- cbind(inData,unlist(gene_z[as.character(inData$gene)]))
  colnames(inData)[3] <- "zsc"
  inData <- inData[match(names(ordered_z), inData$gene),]
  inData$gene <- factor(inData$gene, levels = inData$gene) # orders the factor for plot
  if (useZsc) {
    inData$zsc <- round(inData$zsc,digits = 2)
  } else {
    inData$zsc <- round(inData$zsc,digits = 5)
  }
  
  # makes data to build the legend
  boxPlotLegend <- data.frame(c(26.5,27),c(as.character(inData[nrow(inData)-1,2]),as.character(inData[nrow(inData)-1,2])))
  colnames(boxPlotLegend) <- c("vals","loc")
  inputPlotLegend <- data.frame(c(26.5),c(as.character(inData[nrow(inData)-2,2])))
  colnames(inputPlotLegend) <- c("vals","loc")
  textLegend <- data.frame(c("Serum-Starved","Input Mean"),c(as.character(inData[nrow(inData)-1,2]),as.character(inData[nrow(inData)-2,2])))
  colnames(textLegend) <- c("vals","loc")
  
  # makes the plot
  p <- ggplot(data = inData, aes(x=as.factor(gene), y=expr)) +
    geom_point(color='red', size=3) +
    geom_text(data = inData, aes(x=as.factor(gene), y=20, label=zsc), vjust=0) +
    geom_boxplot(data = distrib, aes(x=as.factor(gene), y=expr), width=.3) +
    geom_point(data = inData, aes(x=as.factor(gene), y=expr), color='red', size=3) +
    coord_flip() +
    geom_line(data = boxPlotLegend, aes(x=as.factor(loc), y=vals), color='black') +
    geom_point(data = inputPlotLegend, aes(x=as.factor(loc), y=vals), color='red', size=3) +
    geom_text(data = textLegend, aes(x=as.factor(loc), y=24, label=vals)) +
    labs(title=paste(NR, " Target Expression", sep=""), 
         subtitle = "Z-Score", x="Target Gene", y="fRMA Expression") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = .7)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  return(p)
}


# This is the main function
CalcEnrich <- function(testSamples){
  test_sample <- testSamples
  
  # keeps track of progress
  num_completed <- 0
  
  # compute scores for all probes (test sample versus serum-starved prior distributions)
  if (ncol(test_sample) == 1) {
    scores <- (test_sample-av)/std
    useZsc <- TRUE
  } else {
    scores <- getDiffex(data,test_sample)
    useZsc <- FALSE
  }

  # creates df to store enrichment p-values
  finPvals <- data.frame(matrix(ncol = length(NRs), nrow = 1))
  colnames(finPvals) <- NRs
  
  finNR <- list() # holds scores for target genes of each NR
  
  for (d in 1:length(NRs)) { # for each nuclear receptor...
    
    # gets a list of target probes for the specified NR
    if (NRs[d] %in% names(allTargets)) {
      target_probes <- allTargets[[NRs[d]]]
    } else {
      next # if NR had less than 15 target genes go to next NR
    }
    
    # select just the probes that are targets of the current NR being evaluated
    test_sample_nr <- scores[target_probes,,drop=F]
    
    # select all other probes that are not listed as targets of the NR
    test_sample_rest <- scores[!row.names(scores)%in%target_probes,,drop=F]
    
    # gets differentially expressed genes from probe z scores
    nrGeneLevelZ <- DiffProbes2Genes(test_sample_nr,useZsc)
    nrGeneLevelZ_mappings <- nrGeneLevelZ[[2]]
    nrGeneLevelZ <- nrGeneLevelZ[[1]]
    
    # generate boxplots of the data
    outPlot <- MakeBoxPlots(test_sample,nrGeneLevelZ,nrGeneLevelZ_mappings,NRs[d],useZsc)
    
    restGeneLevelZ <- DiffProbes2Genes(test_sample_rest,useZsc)
    restGeneLevelZ_mappings <- restGeneLevelZ[[2]]
    restGeneLevelZ <- restGeneLevelZ[[1]]
    
    # # generate boxplots of the data
    # outPlot <- MakeBoxPlots(test_sample,nrGeneLevelZ,nrGeneLevelZ_mappings,NRs[d])
    
    # figure out how many target and non-targets are dysregulated
    nrAboveThresh <- GetNumAboveThreshold(nrGeneLevelZ,useZsc)
    restAboveThresh <- GetNumAboveThreshold(restGeneLevelZ,useZsc)
    
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
    
    # calculates percentage complete
    num_completed <- num_completed + 1
    percent <- round(100*num_completed/15)
    
    incProgress(amount = 1/15,
                detail = paste(percent,"% complete",""))
    
  }
  
  # removes NA values
  finPvals <- finPvals[,!apply(is.na(finPvals), 2, any)]
  
  # appends adjusted p-values
  padj <- p.adjust(finPvals[1,],method = "fdr")
  finPvals <- rbind(finPvals,padj)
  
  # transposes and sorts dataframe
  finPvals <- t(finPvals) 
  finPvals <- finPvals[order(finPvals[,2]),]
  finNR <- as.list(finNR[order(finPvals[,2])]) # orders the list of plots
  # makes rownames into first column
  finPvals <- cbind(Row.Names = rownames(finPvals), finPvals)
  colnames(finPvals) <- c("NR","raw p-value","adj p-value")
  
  # adds gene symbol annotations to differentially expressed probes list
  syms <- as.character(glist[rownames(scores)])
  syms <- replace(syms, syms=="NULL", "")
  scores <- cbind(scores,syms)
  # sorts differentially expressed probes list
  scores <- scores[order(abs(unlist(scores[,1])), decreasing=FALSE),]
  
  # makes probe ID to first column and names columns (for both downloadable tables)
  scores <- cbind(Row.Names = rownames(scores), scores)
  colnames(scores) <- c("Probeset ID","Score","Gene Symbol")
  testSamples <- cbind(Row.Names = rownames(testSamples), testSamples)
  colnames(testSamples)[1] <- c("Probeset ID")
  
  
  return(list(finNR,finPvals,scores,testSamples))
}







