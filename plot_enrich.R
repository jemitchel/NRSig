
library(ggplot2)

results <- readRDS("C:\\Users\\jonat\\Documents\\R\\NRSig-app\\data\\E2_results_eqv.rds")

findb <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(findb) <- c("sample","NR","pval")

for (i in 1:length(results)) {
  pvals <- data.frame(results[[i]][[2]])
  pvals[,3] <- as.numeric(as.character(pvals[,3]))
  for (j in 1:nrow(pvals)) {
    r <- ((i-1)*15)+j
    findb[r,1] <- names(results[i])
    findb[r,2] <- as.character(pvals[j,1])
    
    if (pvals[j,3] < 0.05) {
      findb[r,3] <- -2
    } else if (pvals[j,3] < 0.1) {
      findb[r,3] <- 2
    } else {
      findb[r,3] <- 0
    }
  }
}

NRordered <- list()
for (i in 1:nrow(pvals)) {
  NRordered[[as.character(pvals[i,1])]] <- 0
}
for (i in 1:nrow(findb)) {
  NR <- findb[i,2]
  if (findb[i,3] != 0) {
    if (NR %in% names(NRordered)) {
      NRordered[[NR]] <- NRordered[[NR]] + 1
    }
  }
}

NRordered <- names(NRordered)[order(unlist(NRordered), decreasing = TRUE)]

findb$NR <- factor(findb$NR, levels = NRordered)
findb$sample <- factor(findb$sample, levels = rev(c("ex1","ex2","ex3","ex4","ex5","ex6","ex7","ex8",
                                                "ex9","ex10","ex11","ex12","ex13","ex14","ex15",
                                                "ex16","ex17","ex18","ex19","ex20","ex21")))

colors <- colorRampPalette(c("lightblue", "blue", "gray"))(n=3)
p <- ggplot(data = findb, aes(x = NR, y = sample, fill = pval)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient2(low = colors[2], mid = colors[3], high = colors[1]) +
  scale_x_discrete(position = "top") +
  theme(legend.position="none") +
  coord_equal() + 
  theme(axis.text.x  = element_text(angle=90, vjust=.5, size=10)) +
  theme(axis.text.y  = element_text(angle=0, vjust=0.5, size=6))
p
