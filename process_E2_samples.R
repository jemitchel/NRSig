library(affy)

# gets list of folders each containing different group of replicates to test
# expmnts <- list.dirs("C:/Users/jonat/Documents/R/NRSig-app/data/xeno_E2_experiments")
expmnts <- list.celfiles("C:/Users/jonat/Documents/R/NRSig-app/data/FBS_indv",full.names=TRUE)
# expmnts <- list.celfiles("C:/Users/jonat/Documents/R/NRSig-app/data/xeno_indv",full.names=TRUE)
# expmnts <- list.celfiles("C:/Users/jonat/Documents/R/NRSig-app/data/E2_indv",full.names=TRUE)
# expmnts <- list.files("C:/Users/jonat/Documents/R/NRSig-app/data/E2_indv",full.names=TRUE)



# go to the file with main code to run
source("C:/Users/jonat/Documents/R/NRSig-app/pipeline.R")

# container for results
total_results <- list()

# for (i in 2:length(expmnts)) {
for (i in 1:length(expmnts)) {
# for (i in 37:37) {
# for (i in 4:4) {
  # print(paste("experiment number: ",i-1,""))
  print(paste("experiment number: ",i,""))
  print(expmnts[i])
  
  # get cel files for the replicates
  # files <- list.celfiles(expmnts[i],full.names=TRUE)
  files <- expmnts[i]

  # calculate results
  results <- pipeline(files,NULL)
  # results <- pipeline(files,
  # "C:/Users/jonat/Documents/R/NRSig-app/data/crosshybrid.csv")
  
  # store the results
  total_results[[basename(expmnts[i])]] <- results
}

saveRDS(total_results,"C:/Users/jonat/Documents/R/NRSig-app/data/enr_results/FBS_indv.rds")
