# 
# pipeline <- function(files) {
#   source("C:/Users/jonat/Documents/R/NRSig-app/preprocess.R")
#   samples_matrix <- pre_proc(files)
# 
#   # Computes enrichment
#   source("C:/Users/jonat/Documents/R/NRSig-app/compute_enriched_NRs.R")
#   results <- CalcEnrich(samples_matrix)
#   return(results)
# }
# 
# # cool <- pipeline("C:\\Users\\jonat\\Documents\\Research\\app test folder\\c2_(HG-U133_Plus_2).CEL")
# # c1 <- cool[[1]]
# # c2 <- cool[[2]]
# 
# files <- "C:\\Users\\jonat\\Documents\\Research\\app test folder\\c2_(HG-U133_Plus_2).CEL"
# 
# source("C:/Users/jonat/Documents/R/NRSig-app/preprocess.R")
# samples_matrix <- pre_proc(files)

# Computes enrichment
source("C:/Users/jonat/Documents/R/NRSig-app/compute_enriched_NRs.R")
results <- CalcEnrich(samples_matrix)

