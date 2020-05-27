setwd("E:/Thesis/data/data_camparisions")
#install.packages("readxl")
library(readxl)
file <- "toc_and_elements_all.xlsx"
#install.packages("devtools")
#install.packages("factoextra")
#install.packages("Hmisc")
library("Hmisc")
library("devtools")
library("factoextra")
library(corrplot)
library(rstatix)
library(ggplot2)
library(readxl)
library(gridExtra)
library(ggpubr)
library(ggpmisc)

data <- read_excel(file, sheet = "Sheet1")
data2 <- data[,4:32]
data2 <- data2[,-7]
data2 <- data2[,-16]
data2 <- data2[,-16]
data2 <- data2[,-23]
data2 <- data2[,-23]
df[data2] <- lapply(data2, as.numeric)

round_df(data2, 10)

data2 <- as.matrix(data2)
res2 <- rcorr(as.matrix(data2))
mydata.coeff = res2$r
mydata.p = res2$P

data2 <- data2[,-2]
data2 <- data2[,-4]
data2[is.na(data2)] <- 0.00000000001
 cor.mat <- data2 %>% cor_mat()
pval <- cor.mat %>% cor_get_pval()
cor.mat %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)

test <- cor(as.matrix(data2))
corrplot(test, method = "circle")

write.csv(mydata.coeff,"cor_matrix2.csv")
write.csv(mydata.p,"pval2.csv")
