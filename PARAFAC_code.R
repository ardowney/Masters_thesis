library(dplyr)
library(multiway)
library(staRdom)
library(eemR)
library(parallel)
library(devtools)
library(tidyverse)
library(stats)
library(pracma)
library(tidyr)
library(gtools)
library(tibble)
library(stringr)
library(grDevices)
library(plotly)
library(GGally)

eem_list <- eem_read("percent_corrected", import_function = eem_csv)  #change this line to the desired file path 
dim_min <- 2 # minimum number of components
dim_max <- 7 # maximum number of components
nstart <- 45 # random starts for PARAFAC analysis, models built simulanuously, best selected
cores <- parallel::detectCores(logical=FALSE) # use all cores but do not use all threads
maxit = 7000
ctol <- 10^-7 # tolerance for parafac

#runnigng the model


pfres_comps <- eem_parafac(eem_list, comps = seq(dim_min, dim_max),
                           normalise = TRUE, maxit = maxit, nstart = nstart, ctol = ctol, cores = cores, output = "all")

#Six models are generated within the "pfres_comps" list 
#the following code shoould be ran six times
#each time changing line 34 to analyze one of the six generated models, like so:
#pfmodel <- pfres_comps[[1]]
#pfmodel <- pfres_comps[[2]]
#pfmodel <- pfres_comps[[3]]
#pfmodel <- pfres_comps[[4]]
#pfmodel <- pfres_comps[[5]]
#pfmodel <- pfres_comps[[6]]

pfmodel <- pfres_comps[[1]]

eempf_comp_mat(pfmodel)

eempf_leverage(pfmodel)

leverage <- eempf_leverage(pfmodel)

lev_data <- eempf_leverage_data(leverage) 

pfmodel <- norm2A(pfmodel)

residuals <- eempf_residuals(pfmodel,eem_list)

A_missing(eem_list,pfmodel)

eempf_eemqual(eem_list,pfmodel)

eempf_varimp(pfmodel,eem_list)

comps <- eempf_excomp(pfmodel)


comps <- eempf_excomp(pfmodel,c(1,3))
comps2 <- eempf_excomp(pfmodel,c(4,6))
comps3 <- eempf_bindxc(list(comps, comps2))

factor_table<- eempf_export(pfmodel)

pfres <- pfres_comps[[1]]

eempf_compare(pfmodel)
eempf_leverage_plot(leverage)
#"outliers_1_2"<- eempf_leverage_ident(leverage)
eempf_comp_load_plot(pfmodel)
eempf_comps3D(pfmodel)
eempf_corplot(pfmodel)
eempf_residuals_plot(pfmodel,eem_list)
