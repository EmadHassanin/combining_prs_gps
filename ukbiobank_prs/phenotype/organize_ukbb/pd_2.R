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
               , ggpubr
               , eeptools
               , lubridate)     # add themes and visualization prop. to ggplot   

# 2.0 Importing Files ----

print("source change_me.r file")
source(here::here("phenotype", "organize_ukbb", "change_me.R"))
source(here::here("phenotype", "organize_ukbb", "test.R"))


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
                                              matches("41272"),
                                              matches("41282"),
                                              matches("42032")))

feather_file <- file.path(export_path,"ukb_df","self_reported.feather")
self_reported<- arrow::read_feather(feather_file) %>% as.data.frame()

feather_file <- file.path(export_path,"ukb_df","opcs4.feather")
opcs4 <- arrow::read_feather(feather_file) %>% as.data.frame()

feather_file <- file.path(export_path,"ukb_df","icd10.feather")
icd10 <- arrow::read_feather(feather_file) %>% as.data.frame()

feather_file <- file.path(export_path,"ukb_df","icd9.feather")
icd9 <- arrow::read_feather(feather_file) %>% as.data.frame()

# 3.0 extracting pheno ----

park_icd <- icd10 %>% filter(str_detect(code_icd10 , "G20")) %>% 
   full_join(self_reported %>% filter(code_selfreported == 1262 )) %>%
   mutate(age_selfreported = if_else(age_selfreported < 0 , mean(age_selfreported), age_selfreported )) %>% 
   mutate(age = case_when(
      !is.na(code_icd10) ~ age_icd10,
      !is.na(code_selfreported) ~ age_selfreported,
      !is.na(code_icd10) & !is.na(age_selfreported) ~min(c(age_selfreported,age_icd10))
   )) %>% 
   left_join(ukbb_df %>% select(eid,date_of_parkinsons_disease_report_f42032_0_0)) %>% 
   mutate( age_2 =  time_length(difftime(date_of_parkinsons_disease_report_f42032_0_0,
                                         date_of_birth), "years") ) %>% 
   group_by(eid) %>% 
   slice(which.min(age)) %>% 
   ungroup() %>% 
   select(eid,code_icd10,code_selfreported,age_2,age) %>% 
   mutate(pheno = 1) 

park_icd_final <- 
   ukbb_df %>% select(
      eid,sex_f31_0_0, ethnic_background_f21000_0_0,
      matches("f22009"), sex_chromosome_aneuploidy_f22019_0_0,
      outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
      genetic_sex_f22001_0_0,
      age_when_attended_assessment_centre_f21003_0_0) %>% 
   inner_join(park_icd) %>% 
   filter(ethnic_background_f21000_0_0 == "British") %>% 
   filter(outliers_for_heterozygosity_or_missing_rate_f22027_0_0 != "Yes" | 
             is.na(outliers_for_heterozygosity_or_missing_rate_f22027_0_0)) %>% 
   filter(sex_chromosome_aneuploidy_f22019_0_0 != "Yes" | 
             is.na(sex_chromosome_aneuploidy_f22019_0_0)) %>%
   mutate( sex_check= case_when(
      sex_f31_0_0 == genetic_sex_f22001_0_0 ~ TRUE,
      TRUE ~ FALSE
   )) %>%
   filter(sex_check == TRUE) %>% 
   select(-sex_check) %>% 
   mutate( prev_inc = case_when(
      age_when_attended_assessment_centre_f21003_0_0 > age_2  ~ "prevalent",
      age_when_attended_assessment_centre_f21003_0_0 <= age_2  ~ "incident",
      TRUE  ~ "unknown"
   )) %>% 
   select(
      eid,  age, ethnic_background_f21000_0_0, sex_f31_0_0,
      pheno, matches("f22009"), prev_inc, 
      age_when_attended_assessment_centre_f21003_0_0) 


