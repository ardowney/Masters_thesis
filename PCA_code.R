setwd("E:/Thesis/data/data_camparisions")
file <- "XRF_all_data_molality.xlsx"
library(readxl)
# install.packages("devtools")
# install.packages("factoextra")
library("devtools")
library("factoextra")
data <- read_excel(file, sheet = "molality")


data <- data[,4:28]
data <- data[,2:25]
data <- data[,1:22]
data <- data[,-7]
data <- data[,-4]
data <- data[,-14]
data[is.na(data)] <- 0.00000000001
names(data[, sapply(data, function(v) var(v, na.rm=TRUE)==0)])
res.pca <- prcomp(data, scale = TRUE)
#data <- data[,-5]
#data <- data[,-23]
#data <- data[,-27]
#data2 <- data.frame(t(data[-1]))
#colnames(data2) <- data[,1]
names(data[, sapply(data, function(v) var(v, na.rm=TRUE)==0)])
res.pca <- prcomp(data, scale = TRUE)
fviz_eig(res.pca)
  
#fviz_pca_ind(res.pca,
 #            col.ind = "cos2", # Color by the quality of representation
  #           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
   #          repel = TRUE     # Avoid text overlapping
#)
  fviz_pca_var(res.pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
 # fviz_pca_biplot(res.pca, repel = TRUE,
  #                col.var = "#2E9FDF", # Variables color
   #               col.ind = "#696969"  # Individuals color
  #)
  
  library(factoextra)
  # Eigenvalues
  eig.val <- get_eigenvalue(res.pca)
  eig.val
  
  # Results for Variables
  res.var <- get_pca_var(res.pca)
  res.var$coord          # Coordinates
  res.var$contrib        # Contributions to the PCs
  res.var$cos2           # Quality of representation 
  # Results for individuals
  res.ind <- get_pca_ind(res.pca)
  res.ind$coord          # Coordinates
  res.ind$contrib        # Contributions to the PCs
  res.ind$cos2           # Quality of representation
  

  
  # Helper function 
  #::::::::::::::::::::::::::::::::::::::::
  var_coord_func <- function(loadings, comp.sdev){
    loadings*comp.sdev}
  # Compute Coordinates
  #::::::::::::::::::::::::::::::::::::::::
  loadings <- res.pca$rotation
  sdev <- res.pca$sdev
  var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
  head(var.coord[, 1:4])
  
  var.cos2 <- var.coord^2
  head(var.cos2[, 1:4])
  
  comp.cos2 <- apply(var.cos2, 2, sums
  contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
  var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))

  coord <- as.data.frame(var.coord)
  
  write.csv(eig.val,"eig.csv")
  write.csv(var.coord,"PCA.csv")
  
  library(ggplot2)
  # Basic scatter plot
ggplot(coord, aes(x=PC1, y=PC2)) + geom_point() +
  geom_point() + 
  geom_text(label=rownames(coord), nudge_x=0.02, nudge_y = 0.03)
+
  stat_ellipse(type = "euclid")