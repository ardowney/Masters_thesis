setwd("/Volumes/SAMSUNG_T5/Thesis/data/particle_size/sample_organized")
library(readr)
TAL_16 <- read_csv("TAL-16/TAL-16.csv")
bins <- TAL_16
size <- TAL_16[,1]
res <- apply(array("x", c(nrow("x"), 3, 3)), 3, colMeans = TRUE, na.rm = TRUE)
x<-(TAL_16[,2:188])
y <-rowMeans(x[,1:3])

n <- 1:ncol(x)
unused_df <- matrix(c(n, rep(NA, 3 - ncol(x)%%3)), byrow=TRUE, ncol=3)
unused_df <- data.frame(t(na.omit(unused_df)))
averaged_data <- do.call(cbind, lapply(unused_df, function(i) rowMeans(bins[, i])))

replication_list <- substr(replication_list, 1, nchar(replication_list)-2)
a <- replication_list
b <- a[seq(1, length(a), 3)]
c <- b[40:63]
c <- substr(c, 1, nchar(c)-2)
sample_list <- c(b[1:6],b[8:39],c)
colnames(averaged_data) <- sample_list
rownames(averaged_data) <- size


write.csv(averaged_data, file = "TAL-16_final")
