library(affy)
library(frma)

# samples <- c(list.files("C:/Users/jonat/Documents/R/NRSig-app/data/ss_low_prolif2", full.names=TRUE))

# folders <- list.dirs("C:/Users/jonat/Documents/R/NRSig-app/data/xeno_E2_experiments")
# folders <- list.dirs("C:/Users/jonat/Documents/R/NRSig-app/data/E2_experiments")
folders <- list.dirs("C:/Users/jonat/Documents/R/NRSig-app/data/FBS_experiments")
samples <- c()
for (i in 2:length(folders)) {
  lst <- list.files(folders[[i]], full.names=TRUE)
  samples <- c(samples,lst)
}

batch <- ReadAffy(filenames = samples)
ssData <- frma(batch)
qcreport <- GNUSE(ssData,type = "stats")
data <- exprs(ssData)
print(qcreport)
print(length(samples))

# saveRDS(exprs(ssData),"C:/Users/jonat/Documents/R/NRSig-app/data/serum-starved_lp2.rds")

# res <- barcode(ssData,output = "z-score")
# res["212020_s_at",]
# res["201202_at",]


# thresh1 <- quantile(data2["212020_s_at",])[2]
# thresh2 <- quantile(data2["201202_at",])[2]
# count <- 0
# for (i in 1:ncol(data)) {
#   if (data["212020_s_at",i] < thresh1 || data["201202_at",i] < thresh2) {
#     count <- count + 1
#   } else {
#     print(colnames(data)[i])
#   }
# }
# print(count)

# data <- as.matrix(data)
# hist(data["201202_at",],col=rgb(1,1,0,0.7),freq=F)
# hist(data2["201202_at",],col=rgb(0,1,1,0.4),add=T,freq=F)
# hist(data["212020_s_at",],col=rgb(1,1,0,0.7),freq=TRUE)
# hist(data2["212020_s_at",],col=rgb(0,1,1,0.4),add=T)
# 
# 
# 
# plot.multi.dens <- function(s)
# {
#   junk.x = NULL
#   junk.y = NULL
#   for(i in 1:length(s)) {
#     junk.x = c(junk.x, density(s[[i]])$x)
#     junk.y = c(junk.y, density(s[[i]])$y)
#   }
#   xr <- range(junk.x)
#   yr <- range(junk.y)
#   plot(density(s[[1]]), xlim = xr, ylim = yr, main = "")
#   for(i in 1:length(s)) {
#     lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
#   }
# }
# 
# # the input of the following function MUST be a numeric list
# plot.multi.dens( list(data["201202_at",], data2["201202_at",]))
# plot.multi.dens( list(data["212020_s_at",], data2["212020_s_at",]))

