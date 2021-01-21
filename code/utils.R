# final observation date
last.obs <- as.Date("2/01/2017",format="%m/%d/%Y") 

# Load libraries

library(tidyverse)
library(rlang)
library(magrittr)

library(png)
library(ggalluvial)
library(grid)
library(gridExtra)
library(lattice)
library(RColorBrewer)
library(cowplot)
library(lubridate)


library(survival)
library(survminer)

library(flextable)


# Source codes with functions
source('tidy_data.R')
source('make_newvar.R')
source('make_survanal_var.R')
source('change_time.R')
source('importTidy.R')


source('coverageCalculator.R')
source('create_delayDF.R')
source('KMestim.R')

source("seasonal_func.R")
source("markov_final_VW.R")



# final figure functions
source('fig1_barchart.R')
source('fig2_riverplot.R')
# source('fig3_coverage_timeliness.R')
# source('fig4_survanal.R') 
# source('fig5_late_lollipop.R')
# source('fig6_immunity_gap.R')
# source('fig7_vacc_week.R')
source('fig_vaccinationweek.R')
source('code/raw_code/exploreKM/1_function_getweights.R')


