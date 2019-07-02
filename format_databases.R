library(jsonlite)
library("hgu133plus2.db")
library(annotate)

allTargets <- fromJSON("C:/Users/jonat/Documents/R/NRSig-app/data/TTRUST.json", flatten=TRUE) #use this one for testing
setwd("D:\\Research\\TF Targets")
allTargets2 <- fromJSON("TTRUST_json.txt", flatten=TRUE)

# gets gene symbol annotations for all probes
x <- hgu133plus2SYMBOL
allsymbols <- getSYMBOL(keys(x),"hgu133plus2.db") #list of symbols including probes that don't map

NRs <- c("ESR1", "AR", "NR3C1", "NR3C2", "PGR",
         "NR4A1", "NR4A2", "NR4A3", "NR5A1", "NR5A2", "NR6A1",
         "NR0B1", "NR0B2", "THRA", "THRB", "RARA", "RARB", "RARG",
         "PPARA", "PPARD", "PPARG", "NR1D1", "NR1D2", "RORA", "RORB",
         "RORC", "NR1H4", "NR1H3", "NR1H2", "VDR", "NR1I2",
         "NR1I3", "HNF4A", "HNF4G", "RXRA", "RXRB", "RXRG", "NR2C1",
         "NR2C2", "NR2E1", "NR2E3", "NR2F1", "NR2F2", "NR2F6",
         "ESR2", "ESRRA", "ESRRB", "ESRRG")




# generates list of all gene symbols in 133 plus 2.0 chip
glist <- c()
for (i in 1:length(allsymbols)) {
  g <- allsymbols[[i]]
  if (g %in% glist) {
    # do nothing
  } else {
    glist <- c(glist,g)
  }
}

# eliminates gene symbols in tf-target source not in the chip
for (i in 1:length(NRs)) {
  NR <- NRs[i]
  rems <- c()
  if (NR %in% names(allTargets)) {
    tlist <- allTargets[[NR]]
    if (length(tlist) >= 15) {
      for (j in 1:length(tlist)) {
        gn <- tlist[j]
        if (!(gn %in% glist)) {
          rems <- c(rems,gn)
        }
      }
      
      if (length(rems) > 0) {
        for (j in 1:length(rems)) {
          tlist <- tlist[tlist != rems[j]]
        }
      }
      
      allTargets[[NR]] <- tlist
    } else {
      allTargets <- allTargets[names(allTargets) != NR]
    }
  }
}



all_target_probes <- list()
for (d in 1:length(NRs)) { #for each nuclear receptor...
  # gets a list of target probes for the specified NR
  target_probes <- c()
  target_syms <- allTargets[[NRs[d]]] #first gets list of target genes

  for (i in 1:length(allsymbols)) { #for all probes...
    if (allsymbols[[i]] %in% target_syms) { #if probe maps and maps to a gene in target list
      target_probes <- c(target_probes,names(allsymbols)[i]) #append the probe to a list
    }
  } # CAN REPLACE allTargets with a new json data structure of target probes
  all_target_probes[[NRs[d]]] <- target_probes
}

all_target_probes <- toJSON(all_target_probes)
setwd("C:/Users/jonat/Documents/R/NRSig-app/data")
# write(all_target_probes,"NR_target_probes.json")




