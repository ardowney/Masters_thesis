setwd("F:/Thesis/data/Geology_downey_FS_data")


library("eemR")

library("ggplot2")
library("staRdom")
library("dbplyr")
library("dplyr")

eem_list <- eem_read_csv("CSV_02.21.20/matrix_form")
eem_overview_plot(eem_list, spp=8)
#eem_list <- eem_extend2largest(eem_list, interpolation = 1, extend = FALSE)
#eem_list <- eem_remove_blank(eem_list)
#eem_overview_plot(eem_list, spp=8)
#eem_list <- eem_raman_normalisation2(eem_list, blank = "blank")
#eem_overview_plot(eem_list, spp=8)
#eem_list <- eem_extract(eem_list, c("nano", "miliq", "milliq", "mq", "blank"),ignore_case = TRUE)
remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
remove_scatter_width <- c(15,15,15,15)
eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter, remove_scatter_width = remove_scatter_width)
#eem_overview_plot(eem_list, spp=8)
eem_list <- eem_interp(eem_list, type = 1, extend = FALSE)
eem_overview_plot(eem_list, spp=6)
#test <- c(eem_list[[$ex]])


list_length <- 36
sample <- c()

for (number in 1:list_length) {
  sample <- c(sample, eem_list[[number]][["sample"]])
}
     
extract <- function(EEM){
  y <- EEM[["x"]]
  row <- EEM[["em"]]
  col <- EEM[["ex"]]
  colnames(y) <- (col)
  rownames(y) <- (row)
  return(y)
}

corrected_eem<-lapply(eem_list, extract)

names(corrected_eem) <- sample

setwd("/Volumes/SAMSUNG_T5/Thesis/data/Geology_downey_FS_data/CSV_02.21.20/percent_corrected")
for (file in sample) {
  matrix <- corrected_eem[[file]]
  write.csv(matrix, file = paste0(file, "corrected.csv"))
}
  
files <- list.files()


for (file_name in files){
  x_data <- read.csv(file_name)
  rows <- x_data[,1]
  col <- colnames(x_data)
  col <- sub("X", "", col)
  colnames(x_data) <- (col)
  row.names(x_data) <- rows
  x_data <- x_data[,-1]
  max <- max(x_data)
  x_data2 <- (x_data)/(max)
  assign(file_name, x_data2)
  write.csv(x_data2, file = sub(".csv","",paste0(file_name, "percent.csv")))
}


setwd("/Volumes/SAMSUNG_T5/Thesis/data/Geology_downey_FS_data/CSV_02.21.20")
eem_list <- eem_read_csv("percent_corrected")
eem_overview_plot(eem_list, spp=8)
#eem_list <- eem_extend2largest(eem_list, interpolation = 1, extend = FALSE)
#eem_list <- eem_remove_blank(eem_list)
#eem_overview_plot(eem_list, spp=8)
#eem_list <- eem_raman_normalisation2(eem_list, blank = "blank")
#eem_overview_plot(eem_list, spp=8)
#eem_list <- eem_extract(eem_list, c("nano", "miliq", "milliq", "mq", "blank"),ignore_case = TRUE)
remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
remove_scatter_width <- c(15,15,15,15)
eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter, remove_scatter_width = remove_scatter_width)
#eem_overview_plot(eem_list, spp=8)
eem_list <- eem_interp(eem_list, type = 1, extend = FALSE)
eem_overview_plot(eem_list, spp=6)
#test <- c(eem_list[[$ex]])


peaks <- eem_coble_peaks(eem_list)

FI <- eem_fluorescence_index(eem_list, verbose = TRUE)
HI <- eem_humification_index(eem_list, scale = FALSE, verbose = TRUE)
BI <- eem_biological_index(eem_list)


setwd("/Volumes/SAMSUNG_T5/Thesis/data/Geology_downey_FS_data/CSV_02.21.20/summary_parameters")

write.csv(peaks, "peaks.csv")
write.csv(BI, "BI.csv")
write.csv(FI, "FI.csv")
write.csv(HI, "HI.csv")
