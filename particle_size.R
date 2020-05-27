setwd("/Volumes/SAMSUNG_T5/Thesis/data/particle_size/sample_organized")
size_file <- "particle_size_bins.txt"
size <- read.delim(size_file, header = FALSE)
setwd("/Volumes/SAMSUNG_T5/Thesis/data/particle_size/sample_organized/TAL-15")
files <- list.files()
bins <- data.frame(size)
row.names(bins) <- bins[1:92,]
bins <- bins[1:92,0]
for (xfile_name in files) {
  x_data <- read.delim(paste0(xfile_name), header = TRUE)
  assign(xfile_name, x_data) 
  bins <- cbind(bins, get(xfile_name))
}

colnames(bins) <- files
write.csv(bins, file = "TAL-15.csv")
