setwd("F:/Thesis/data/Geology_downey_FS_data")
files <- list.files()


for (file_name in files){
  setwd("F:/Thesis/EEM/3-D_spectra/sed_extracts_101119/test")
  x_data <- read.csv(paste(file_name))
max <- max(x_data)
x_data2 <- (x_data)/(max)
assign(file_name, x_data2)
setwd("F:/Thesis/EEM/3-D_spectra/sed_extracts_101119/percent")
write.csv(x_data, file = sub(".csv","",paste0(file_name, "percent.csv")))
}
col <- colnames(x_data)
col <- sub("X", "", col)
colnames(x_data) <- (col)
