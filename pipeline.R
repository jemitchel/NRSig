source("C:/Users/jonat/Documents/R/NRSig-app/preprocess.R")

# A list of files is collected for analysis from GUI (user selected)
flist <- c("D:\\Research\\test folder\\c3_(HG-U133_Plus_2).CEL","D:\\Research\\test folder\\c2_(HG-U133_Plus_2).CEL")
samples_matrix <- pre_proc(flist)

# Compute enrichment
source("C:/Users/jonat/Documents/R/NRSig-app/compute_enriched_NRs.R")
CalcEnrich(samples_matrix)
