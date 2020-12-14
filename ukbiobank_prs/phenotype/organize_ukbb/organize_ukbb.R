# Organize uk_biobank phenotype data ----

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
             #  , arrow       # apache arrow provides powerful package for sharing data
               , data.table  # exporting data as .txt
               , ggpubr      # add themes and visualization prop. to ggplot
               , doParallel
               , here
               ,glmnet)    

#install.packages("arrow", repos = c("https://dl.bintray.com/ursalabs/arrow-r"))
                 
#library("arrow")
# 2.0 Importing Files ----

print("source change_me.r file")
source(here::here("phenotype","organize_ukbb","change_me.R"))
source(here::here("phenotype","organize_ukbb","test.R"))

feather_file <- file.path(export_path,"ukb_df","ukb_df.rds")


print("time taken to import all datasets")

tic()
ukbb_pc_df <- arrow::read_feather(feather_file,
                                  col_select = c(matches("22009"), col_vars$col_name))

ukbb_df <- arrow::read_feather(feather_file,
                               col_select = c(col_vars$col_name,  matches("age_at_cancer_diagnosis") ,
                                              matches("cancer_code_selfreported"),
                                              matches("noncancer_illness_code_selfreported"),
                                              matches("type_of_cancer"),
                                              matches("interpolated_age_of_participant"),
                                              matches("date_of_cancer_diagnosis"),
                                              matches("^diagnoses.*icd"),
                                              matches("^diagnoses.*icd"),
                                              behaviour_of_cancer_tumour_f40012_0_0,
                                              matches("date_of_attending_assessment_centre"),
                                              matches("22009"), matches("41280"),matches("41262"),
                                              matches("41281"), matches("41263"),
                                              matches("month_of_birth_f52"),
                                              matches("year_of_birth_f34_0_0"),
                                              matches("illnesses_of"),
                                              matches("f87"),
                                              matches("30770")))


toc()

tic()
print("time taken to import all diagnosis based on icd codes")

feather_file <- file.path(export_path,"ukb_df","icdcodes_all.rds")
icdcodes_all <- arrow::read_feather(feather_file) %>% as.data.frame()

toc()

fam_files <- list.files(fam_files_path, full.names = T)
#fam_df <- fread(file.path(fam_files_path,"ukb_chr1.fam" ))

# diagnosis data

feather_file <- file.path(export_path,"ukb_df","self_reported.feather")
self_reported<- arrow::read_feather(feather_file) %>% as.data.frame()

feather_file <- file.path(export_path,"ukb_df","opcs4.feather")
opcs4 <- arrow::read_feather(feather_file) %>% as.data.frame()

feather_file <- file.path(export_path,"ukb_df","icd10.feather")
icd10 <- arrow::read_feather(feather_file) %>% as.data.frame()

feather_file <- file.path(export_path,"ukb_df","icd9.feather")
icd9 <- arrow::read_feather(feather_file) %>% as.data.frame()


# 3.0 Diagnosis ----

tic()
print("time taken to extract phenotypic information to diagnosed samples")



pheno_icd_ukbtools <- 
   inner_join( ukbb_df  %>% 
                  select( col_vars$col_name,  matches("age_at_cancer_diagnosis") ,
                          matches("interpolated_age_of_participant_when_cancer_first"),
                          matches("date_of_cancer_diagnosis"),
                          behaviour_of_cancer_tumour_f40012_0_0,
                          matches("date_of_attending_assessment_centre")),
               icdcodes_all %>% 
                  filter(str_detect(icd9_code , identify_icd$icd9)   | 
                            str_detect(icd10_code ,identify_icd$icd10)),
               by = c( "eid" = "sample" )) %>% 
   full_join(self_reported_cancer_diagnosis(ukbb_df, 1002)) 

toc()

pheno_icd_all<- 
  inner_join( ukbb_df  %>% 
                 select( eid, col_vars$col_name,  matches("age_at_cancer_diagnosis") ,
                         matches("interpolated_age_of_participant_when_cancer_first"),
                         matches("date_of_cancer_diagnosis"),
                         behaviour_of_cancer_tumour_f40012_0_0,
                         matches("date_of_attending_assessment_centre")),
              icdcodes_all %>% 
                filter(str_detect(icd9_code , identify_icd$icd9)   | 
                         str_detect(icd10_code ,identify_icd$icd10)),
              by = c( "eid" = "sample" ))  %>% 
   full_join(register_cancer_diagnosis(ukbb_df,"^(C50[0-9])|^(174[0-9])")) %>% 
   full_join(self_reported_cancer_diagnosis(ukbb_df, 1002))
#   inner_join(fam_df %>% select(V1), by = c( "eid" = "V1" ))

pheno_icd_registered_1 <- 
   
   inner_join( ukbb_df  %>% 
                  select( col_vars$col_name,  matches("age_at_cancer_diagnosis") ,
                          matches("interpolated_age_of_participant_when_cancer_first"),
                          matches("date_of_cancer_diagnosis"),
                          behaviour_of_cancer_tumour_f40012_0_0,
                          matches("date_of_attending_assessment_centre")),
               register_cancer_diagnosis(ukbb_df,"^(C50[0-9])|^(174[0-9])"),
               by = c( "eid" = "eid"  )) 

pheno_icd_registered_2 <- 
   
   inner_join( ukbb_df  %>% 
                  select( col_vars$col_name,  matches("age_at_cancer_diagnosis") ,
                          matches("interpolated_age_of_participant_when_cancer_first"),
                          matches("date_of_cancer_diagnosis"),
                          behaviour_of_cancer_tumour_f40012_0_0,
                          matches("date_of_attending_assessment_centre")),
               self_reported_cancer_diagnosis(ukbb_df, 1002),
               by = c( "eid" = "eid"  ))





pheno_icd_registered <- bind_rows(pheno_icd_registered_1,pheno_icd_registered_2)




breast_fem_white <- pheno_icd_registered %>% 
   mutate( min_age_at_cancer_diagnosis =  apply(pheno_icd_registered %>% select( contains("age_at_cancer_diagnosis"),
                                                                         matches("interpolated_age_of_participant_when_cancer_first")),
                                                1,min, na.rm = TRUE),
           min_age_at_cancer_diagnosis = ifelse(is.infinite(min_age_at_cancer_diagnosis ), 
                                                NA, min_age_at_cancer_diagnosis )) %>%
   mutate( prev_inc = case_when(
      age_when_attended_assessment_centre_f21003_0_0 > min_age_at_cancer_diagnosis  ~ "prevalent",
      age_when_attended_assessment_centre_f21003_0_0 < min_age_at_cancer_diagnosis  ~ "incident",
      age_when_attended_assessment_centre_f21003_0_0 == min_age_at_cancer_diagnosis ~ "unknown"
   )) %>% 
   filter(!is.na(behaviour_of_cancer_tumour_f40012_0_0)) %>% 
   filter(sex_f31_0_0 == "Female" & ethnic_background_f21000_0_0 == "British") %>% 
   filter(outliers_for_heterozygosity_or_missing_rate_f22027_0_0 != "Yes" | 
          is.na(outliers_for_heterozygosity_or_missing_rate_f22027_0_0)) %>% 
   filter(sex_chromosome_aneuploidy_f22019_0_0 != "Yes" | 
             is.na(sex_chromosome_aneuploidy_f22019_0_0)) %>%
   mutate( sex_check= case_when(
      sex_f31_0_0 == genetic_sex_f22001_0_0 ~ TRUE,
      TRUE ~ FALSE
   )) %>%
   filter(sex_check == TRUE) %>% 
   select(-sex_check) 
   

breast_fem_indian <- pheno_icd_registered %>% 
   mutate( min_age_at_cancer_diagnosis =  apply(pheno_icd_registered %>% select( contains("age_at_cancer_diagnosis"),
                                                                                 matches("interpolated_age_of_participant_when_cancer_first")),
                                                1,min, na.rm = TRUE),
           min_age_at_cancer_diagnosis = ifelse(is.infinite(min_age_at_cancer_diagnosis ), 
                                                NA, min_age_at_cancer_diagnosis )) %>%
   mutate( prev_inc = case_when(
      age_when_attended_assessment_centre_f21003_0_0 > min_age_at_cancer_diagnosis  ~ "prevalent",
      age_when_attended_assessment_centre_f21003_0_0 < min_age_at_cancer_diagnosis  ~ "incident",
      age_when_attended_assessment_centre_f21003_0_0 == min_age_at_cancer_diagnosis ~ "unknown"
   )) %>% 
   filter(sex_f31_0_0 == "Female" & ethnic_background_f21000_0_0 == "Indian") %>% 
   filter(outliers_for_heterozygosity_or_missing_rate_f22027_0_0 != "Yes" | 
             is.na(outliers_for_heterozygosity_or_missing_rate_f22027_0_0)) %>% 
   filter(sex_chromosome_aneuploidy_f22019_0_0 != "Yes" | 
             is.na(sex_chromosome_aneuploidy_f22019_0_0)) %>%
   mutate( sex_check= case_when(
      sex_f31_0_0 == genetic_sex_f22001_0_0 ~ TRUE,
      TRUE ~ FALSE
   )) %>%
   filter(sex_check == TRUE) %>% 
   dplyr::select(-sex_check)

#toc()

# 4.0 exporting files  ----

#tic()
#print("time taken to write the data")

#fwrite(pheno_icd, file = file.path(export_path,
#                                    paste(identify_icd$phenotype , '_phenotype.txt', sep ="")), 
#           sep="\t", row.names = FALSE)
#toc()

#tic()
#print("time taken to save the data as .rds")
#saveRDS(pheno_icd, file = file.path(export_path,
#                                    paste(identify_icd$phenotype , '_phenotype.rds', sep ="")))

#toc()











