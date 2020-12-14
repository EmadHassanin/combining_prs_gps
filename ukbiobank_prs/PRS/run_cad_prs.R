# 1.0 Load libraries ----
print("install and loading dependencies packages")

installed.packages("pacman") 
library("pacman")

#devtools::install_github("kenhanscombe/ukbtools", build_vignettes = TRUE, dependencies = TRUE, force = TRUE)

pacman::p_load(tidyverse     # tidy daya
               , ukbtools    # workhorse package for organizing ukbb
               , tictoc      # checking time 
               , furrr       # functional programing with map()
               , future      # parallel computing
               , arrow       # apache arrow provides powerful package for sharing data
               , data.table  # exporting data as .txt
               , ggpubr      # add themes and visualization prop. to ggplot
               ,plotly       # interactive plots
               ,caret        # partition data
               #,rcompanion   # calculating R2
               ,pROC         # calculating and plot AUC 
               ,plotROC      # calculating and plot ggplot
               ,gridExtra    # add table to ggplot
               , here
)    

#import files
#source(here::here("phenotype","organize_ukbb","change_me.R"))
#source(here::here("phenotype","organize_ukbb","test.R"))
#source(here::here("phenotype","organize_ukbb","organize_ukbb.R"))
#source(here::here("phenotype","organize_ukbb","cad.R"))
source(here::here("PRS","cad_prs.R"))
