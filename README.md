# Physical and chemical properties of clastic cave sediment from the northern karst region of Puerto Rico 

## R scripts include:
### Organizing raw particle size data from a Beckman Coulter single wavelength LS13-320 particle size analyzer from the instrument into a csv format

### post-processing fluorescence spectroscopy data from a Horiba Fluorolog-3

### PARAFAC modeling of EEM data

### PCA analysis

### Pearson Correlation Matrix 

## Prerequisites:
### Communication with Unix Shell

##Necessary R packages:
	###Particle size -### install.packages("tidyverse")
      			 ###  install.packages("dplyr")
      			 ###  install.packages("readxl") 
      			 ###  install.packages("ggplot2")
      			 ###  install.packages("colorspace")
      			 ###  install.packages("scales")
## Fluorescence data – ### install.packages("eemR")
		       ### install.packages("ggplot2")
		       ### install.packages("staRdom")
install.packages("dbplyr")
install.packages("dplyr")
###PARAFAC model - install.packages(dplyr)
 install.packages(multiway)
 install.packages(staRdom)
 install.packages(eemR)
 install.packages(parallel)
 install.packages(devtools)
 install.packages(tidyverse)
 install.packages(stats)
 install.packages(pracma)
 install.packages(tidyr)
 install.packages(gtools)
 install.packages(tibble)
 install.packages(stringr)
 install.packages(grDevices)
 install.packages(plotly)
 install.packages(GGally)

###PCA – install.packages(readxl)
 install.packages("devtools")
 install.packages("factoextra")

###Pearson Correltion Matrix – install.packages(readxl)
   install.packages("Hmisc")
   install.packages("devtools")
   install.packages("factoextra")
   install.packages(corrplot)
   install.packages(rstatix)
   install.packages(ggplot2)
   install.packages(readxl)
   install.packages(gridExtra)
   install.packages(ggpubr)
   install.packages(ggpmisc)

##script name and function:

###particle size:
cut_loop.sh – Unix shell script organizing raw data from the instrument. transforms proprietary file extension to flat text file (.txt) and sends to a new organized directory. 

###Particle size.R – applies particle size bins from the instrument to the organized data created by cut_loop.sh. 

###Averaging_script.R – takes the average of the three replicates ran for each subsample. 

###Fluorescence data: 
Scatter_code.R – plots EEM data stored as a csv. Removes 15 nm 1st and 2nd order Rayleigh scattering. Smooths the EEM by interpreting intensity through the Rayleigh regions. Sends corrected files to a modified csv file. 

###Percent_function.R – normalizes corrected EEM to the max intensity of each scan. 

###PARAFAC model:
Model_functions.R – Functions needed to run the PARAFAC model

PARAFAC_code.R – script for creating 2-7 component models using user defined EEM data. 

###PCA: 
PCA.R – creates PCA output for user defined data. Outputs include biplots, variance plots, eigenvalues, and standard PCA scatter plots. 

Pearson Correlation Matrix:
Cor_matrix_function.R – creates correlation matrix with R and p values exported to a csv file. R values and p values plots are included within the code.  
