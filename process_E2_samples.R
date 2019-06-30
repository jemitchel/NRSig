
# gets list of folders each containing different group of replicates to test
# expmnts <- list.dirs("C:/Users/jonat/Documents/R/NRSig-app/data/E2_experiments")
expmnts <- list.dirs("C:/Users/jonat/Documents/R/NRSig-app/data/xeno_E2_experiments")


# go to the file with main code to run
source("C:/Users/jonat/Documents/R/NRSig-app/pipeline.R")

# container for results
total_results <- list()

for (i in 2:length(expmnts)) {
  print(paste("experiment number: ",i-1,""))
  print(expmnts[i])
  
  # get cel files for the replicates
  files <- list.celfiles(expmnts[i],full.names=TRUE) 

  # calculate results
  # results <- pipeline(files,NULL)
  results <- pipeline(files,
                      "C:/Users/jonat/Documents/R/NRSig-app/data/crosshybrid_fives.csv")
  
  # store the results
  total_results[[basename(expmnts[i])]] <- results
}

# saveRDS(total_results,"C:/Users/jonat/Documents/R/NRSig-app/data/xeno_results_qc.rds")
