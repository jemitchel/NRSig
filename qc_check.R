library(affy)
library(frma)

# samples <- c(list.files("C:/Users/jonat/Documents/R/NRSig-app/data/ss_samples", full.names=TRUE))

folders <- list.dirs("C:/Users/jonat/Documents/R/NRSig-app/data/xeno_E2_experiments")
samples <- c()
for (i in 2:length(folders)) {
  lst <- list.files(folders[[i]], full.names=TRUE)
  samples <- c(samples,lst)
}

batch <- ReadAffy(filenames = samples)
ssData <- frma(batch)
qcreport <- GNUSE(ssData,type = "stats")
qcreport 
# print(length(samples))
# 
# saveRDS(exprs(ssData),"C:/Users/jonat/Documents/R/NRSig-app/data/serum-starved_qc.rds")
# 
