# 1.0 Load libraries ----
print("install and loading dependencies packages")


installed.packages("pacman") 

library("pacman")

#devtools::install_github("kenhanscombe/ukbtools", build_vignettes = TRUE, dependencies = TRUE, force = TRUE)

pacman::p_load(tidyverse     # tidy daya
               ,
               ukbtools    # workhorse package for organizing ukbb
               ,
               tictoc      # checking time
               ,
               furrr       # functional programing with map()
               ,
               future      # parallel computing
               ,
               arrow       # apache arrow provides powerful package for sharing data
               ,
               data.table  # exporting data as .txt
               ,
               ggpubr
               , 
               
               eeptools)     # add themes and visualization prop. to ggplot   

# 2.0 Importing Files ----

print("source change_me.r file")
source(here::here("phenotype", "organize_ukbb", "change_me.R"))
source(here::here("phenotype", "organize_ukbb", "test.R"))

feather_file <- file.path(export_path, "ukb_df", "ukb_df.rds")


print("time taken to import all datasets")



ukbb_df_park <- arrow::read_feather(
   feather_file,
   col_select = c(
      col_vars$col_name,
      matches("noncancer_illness_code_selfreported"),
      matches("^diagnoses.*icd"),
      matches("^diagnoses.*icd"),
      matches("date_of_attending_assessment_centre"),
      matches("22009"),
      matches("illnesses_of"),
      matches("20009"),
      matches("birth"),
      matches("41280"),
      
   )
)



feather_file <- file.path(export_path, "ukb_df", "icdcodes_all.rds")
icdcodes_all <- arrow::read_feather(feather_file) %>% as.data.frame()


park_icd_ukbtools_1 <-
   inner_join(
      ukbb_df_park,
      icdcodes_all %>%
         filter(str_detect(icd10_code , "G20")),
      by = c("eid" = "sample")
   )


park_icd_ukbtools_2 <- ukbb_df_park %>% 
   inner_join(self_reported_non_cancer_diagnosis(ukbb_df, 1288))  

park_icd_ukbtools <- bind_rows(park_icd_ukbtools_1,park_icd_ukbtools_2)

tr <- park_icd_ukbtools<- 
   icd10 %>% filter(str_detect(code_icd10 , "G20")) %>% 
   full_join(self_reported %>% filter(code_selfreported == 1288 ))



park_icd_ukbtools %>%
   mutate(
      month_of_birth_f52_0_0 = match(month_of_birth_f52_0_0, month.name),
      date_of_birth = as.Date(paste(year_of_birth_f34_0_0, month_of_birth_f52_0_0, "01", sep = "-")),
      min_date_at_diagnosis =  apply(park_icd_ukbtools %>% select(contains("41280")),
                                     1, min, na.rm = TRUE),
      min_date_at_diagnosis = ifelse(is.infinite(min_date_at_diagnosis),NA,min_date_at_diagnosis),
      age = time_length(difftime(min_date_at_diagnosis, date_of_birth), "years"),
      min_age_at_diagnosis =  apply(park_icd_ukbtools %>% select(contains("20009")),
                                    1, min, na.rm = TRUE),
      min_age_at_diagnosis = ifelse(is.infinite(min_age_at_diagnosis),
                                    NA, min_age_at_diagnosis)) %>%
   mutate(type = case_when(
      (is.na(icd10_code) == FALSE & is.na(code) == FALSE) ~ "both",
      is.na(icd10_code) == FALSE ~ "icd_code10" ,
      is.na(code) == FALSE ~ "self_reported"
   )) %>%
   mutate(pheno = 1) %>%
   
   #distinct(eid, .keep_all = TRUE) %>% 

   #filter(ethnic_background_f21000_0_0 == "British") %>% 
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
      age_when_attended_assessment_centre_f21003_0_0 > age  ~ "prevalent",
      age_when_attended_assessment_centre_f21003_0_0 < age  ~ "incident",
      age_when_attended_assessment_centre_f21003_0_0 == age ~ "unknown"
   )) %>% 
   select(
      eid,  age, pheno, genetic_principal_components_f22009_0_1, ethnic_background_f21000_0_0,
      genetic_principal_components_f22009_0_2, type,prev_inc) -> park_ukbtools

control_ids <- 
   ukbb_df %>% dplyr::select(eid,sex_f31_0_0,ethnic_background_f21000_0_0,
                             outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
                             sex_chromosome_aneuploidy_f22019_0_0,
                             genetic_sex_f22001_0_0,
                             age_when_attended_assessment_centre_f21003_0_0,
                             matches("f22009")) %>%
  # filter(ethnic_background_f21000_0_0 == "British") %>%
   filter(outliers_for_heterozygosity_or_missing_rate_f22027_0_0 != "Yes" |
             is.na(outliers_for_heterozygosity_or_missing_rate_f22027_0_0)) %>%
   filter(sex_chromosome_aneuploidy_f22019_0_0 != "Yes" |
             is.na(sex_chromosome_aneuploidy_f22019_0_0)) %>%
   mutate( sex_check= case_when(
      sex_f31_0_0 == genetic_sex_f22001_0_0 ~ TRUE,
      TRUE ~ FALSE
   )) %>%
   filter(sex_check == TRUE) %>%
   dplyr::select(-sex_check) %>%
   anti_join(park_ukbtools, by = c("eid")) %>% 
   dplyr::select(eid,age_when_attended_assessment_centre_f21003_0_0,ethnic_background_f21000_0_0,
                 matches("f22009"))  %>%
   rename(age = age_when_attended_assessment_centre_f21003_0_0) %>% 
   distinct(.keep_all = TRUE) %>%
   mutate(pheno = 0,
          type = "control",
          prev_inc = "control")

all_ids <- park_ukbtools %>%
   bind_rows(control_ids)

fwrite(all_ids, file = file.path(export_path, "output","ukbb_pheno_test_output","park_ukbb.txt"), 
           sep="\t", row.names = FALSE)

emad <- fread(file = file.path(export_path, "output","ukbb_pheno_test_output","park_ukbb.txt"))

park_icd_ukbtools_1 <-
   inner_join(
      ukbb_df_park,
      icd10_diagnosis %>% select(eid,code,age) %>% 
         filter(str_detect(code , "G20")),
      by = c("eid" = "eid")
   )

park_icd_ukbtools_2 <- inner_join(
   ukbb_df_park,
   self_reported_diagnosis %>% select(eid, code_selfreported, age_selfreported) %>% 
      filter(str_detect(code_selfreported , "1288")),
   by = c("eid" = "eid")
) 


park_icd_ukbtools %>% 
   mutate(type = case_when(
      (is.na(code) == FALSE & is.na(code_selfreported) == FALSE) ~ "both",
      is.na(code) == FALSE ~ "icd_code10" ,
      is.na(code_selfreported) == FALSE ~ "self_reported"
   )) %>%
   mutate(pheno = 1) %>% 
   mutate(age = ifelse( type == "self_reported" , age_selfreported, age))   %>% 
   #distinct(eid, .keep_all = TRUE) %>% 
   
   #filter(ethnic_background_f21000_0_0 == "British") %>% 
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
      age_when_attended_assessment_centre_f21003_0_0 > age  ~ "prevalent",
      age_when_attended_assessment_centre_f21003_0_0 < age  ~ "incident",
      age_when_attended_assessment_centre_f21003_0_0 == age ~ "unknown"
   )) %>% 
   select(
      eid,  age, pheno, matches("f22009"), type,prev_inc) 

