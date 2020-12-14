# Diagnosis of uk_biobank phenotypes according to icd_codes 9 and 10 versions ----

# 1.0 Load libraries ----
print("install and loading dependencies packages")

installed.packages("pacman") 
library("pacman")

#devtools::install_github("kenhanscombe/ukbtools", build_vignettes = TRUE, dependencies = TRUE, force = TRUE)

pacman::p_load(tidyverse, ukbtools, tictoc, furrr, future)

# 2.0 Importing Files ----
print("source change_me.r file")
source(here::here("organize_ukbb","change_me.R"))

tic()
interested_ukb_df <- readRDS(file =  file.path( data_path, file_name)) %>% 
  select(col_vars$col_name, identify_columns, 
         matches("^diagnoses.*icd9"), matches("^diagnoses.*icd10")
  )

print("time taken to import all datasets")
toc()

# 3.0 Diagnosis ----

no_cores <- availableCores() - 1
plan(multisession, workers = no_cores)
tic()

print("time taken to extract diagnosed samples")
options(future.globals.maxSize = 3000 * 1024^2)

icd_diagnosis_data_1 %<-% ukb_icd_diagnosis(interested_ukb_df  , id = interested_ukb_df$eid , icd.version = 9)
icd_diagnosis_data_2 %<-% ukb_icd_diagnosis(interested_ukb_df , id = interested_ukb_df$eid , icd.version = 10)

icd_diagnosis_data <- full_join(icd_diagnosis_data_1,icd_diagnosis_data_2,by = c("sample" = "sample")) %>% 
  mutate(sample = as.integer(sample) ) %>% 
  rename(icd9_code = code.x , icd9_meaning = meaning.x,
         icd10_code = code.y , icd10_meaning = meaning.y)
toc()

# 4.0 Export files ----

tic()
print("time taken to save the data as .rds")
saveRDS(icd_diagnosis_data, file = file.path(export_path,'icdcodes_all.rds'))

toc()

feather_file <- file.path(export_path,"ukb_df","icdcodes_all.rds")
arrow::write_feather(uicd_diagnosis_data, feather_file)