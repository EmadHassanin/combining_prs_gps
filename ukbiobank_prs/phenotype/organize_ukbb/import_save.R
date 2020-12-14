# Import and save uk_biobank phenotype datasets ----
# 1.0 Load libraries ----

print("install and loading dependencies packages")

installed.packages("pacman") 
library("pacman")


pacman::p_load(tidyverse, ukbtools, tictoc, furrr, future)

#devtools::install_github("kenhanscombe/ukbtools", build_vignettes = TRUE, dependencies = TRUE, force = TRUE)

# 2.0 Importing Files ----
print("source change_me.r file")
source(here::here("phenotype" ,"organize_ukbb","change_me.R"))

data_path <- "/home/opc/storage/uk_biobank/phenotype_data"
files <- c("ukb39215","ukb41120")

print("imports datasets in parallel using furrr package")

no_cores <- availableCores() - 1
plan(multisession, workers = no_cores)

tic()
ukbb_df <- files %>% 
  future_map(ukb_df, path = data_path,  n_threads = "max") %>% 
  reduce(right_join)
print("time taken for imporing datat")
toc()


# 3.0 saving Files ----

print("save datasets as .rds")

tic()
saveRDS(ukbb_df, file = file.path(export_path,'ukb_df.rds'))
toc()

print("save datasets as .txt")
tic()
write.table(ukbb_df, file = file.path(export_path,'ukb_df.txt'), 
            sep = "\t" , row.names = FALSE)
toc()

tic()
readRDS(file = "my_data.rds")
feather_file <- file.path(export_path,"ukb_df","ukb_df.rds")
arrow::write_feather(ukbb_df, feather_file)
toc()

print("time taken to import all datasets")

# print("load datasets .rds")
# tic()
# ukbiobank_500 <- readRDS(file = file.path(export_path,'ukb_df_500.rds')) %>% 
#   select( col_vars$col_name, identify_columns,
#          matches("^diagnoses.*icd9"), matches("^diagnoses.*icd10")
#   )
# toc()
# dim(ukbiobank_500)
