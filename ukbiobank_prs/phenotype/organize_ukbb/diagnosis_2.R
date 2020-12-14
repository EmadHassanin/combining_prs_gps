# self reported cases ----

self_reported_codes <- ukbb_df %>% 
   dplyr::select(eid, matches("noncancer_illness_code_selfreported") ) %>% 
   gather("non_cancer_code_selfreported_version" , "code_selfreported", - eid) %>% 
   mutate(new = str_remove(non_cancer_code_selfreported_version,
                           "noncancer_illness_code_selfreported_f20002_")) 

self_reported_age_diagnosis <- ukbb_df %>% 
   dplyr::select(eid, matches("20009") ) %>% 
   gather("age_version" , "age_selfreported", - eid) %>%
   mutate(new = str_remove(age_version,
                           "interpolated_age_of_participant_when_noncancer_illness_first_diagnosed_f20009_"))

self_reported_codes  %>% 
   left_join(self_reported_age_diagnosis ) -> self_reported_diagnosis

self_reported_diagnosis %>% 
   filter(!is.na(code_selfreported)) -> self_reported

feather_file <- file.path(export_path,"ukb_df","self_reported.feather")
arrow::write_feather(self_reported, feather_file)
# ICD 10 + main ICD 10 ----

birth_ukbb_df <- ukbb_df %>% select(eid, month_of_birth_f52_0_0, year_of_birth_f34_0_0) %>%
   mutate(
      month_of_birth_f52_0_0 = match(month_of_birth_f52_0_0, month.name),
      date_of_birth = as.Date( paste(year_of_birth_f34_0_0, month_of_birth_f52_0_0, "01", sep = "-" ))) %>% 
   select(eid,date_of_birth )

icd10_codes <- ukbb_df %>% 
   dplyr::select(eid, matches("41270") ) %>% 
   gather("icd10_diagnosis" , "code_icd10", - eid) %>%  
   mutate(new = str_remove(icd10_diagnosis,
                           "diagnoses_icd10_f41270_")) 

icd10_age_diagnosis <- ukbb_df %>% 
   dplyr::select(eid, matches("41280") ) %>% 
   gather("date_version" , "date", - eid) %>% 
   mutate(new = str_remove(date_version,
                           "date_of_first_inpatient_diagnosis_icd10_f41280_")) %>% 
   left_join(birth_ukbb_df) %>% 
   mutate(age_icd10 = time_length(difftime(date, date_of_birth), "years"))

icd10_codes  %>% 
   left_join(icd10_age_diagnosis) -> icd10_diagnosis

icd10_diagnosis %>% 
   filter(!is.na(code_icd10)) -> icd10

feather_file <- file.path(export_path,"ukb_df","icd10.feather")
arrow::write_feather(icd10, feather_file)
# ICD 9 + main ICD 9 ----

icd9_codes <- ukbb_df %>% 
   dplyr::select(eid, matches("41271") ) %>% 
   gather("icd9_diagnosis" , "code_icd9", - eid) %>%  
   mutate(new = str_remove(icd9_diagnosis,
                           "diagnoses_icd9_f41271_")) 

icd9_age_diagnosis <- ukbb_df %>% 
   dplyr::select(eid, matches("41281") ) %>% 
   gather("date_version" , "date", - eid) %>% 
   mutate(new = str_remove(date_version,
                           "date_of_first_inpatient_diagnosis_icd9_f41281_")) %>% 
   left_join(birth_ukbb_df) %>% 
   mutate(age_icd9 = time_length(difftime(date, date_of_birth), "years"))

icd9_codes  %>% 
   left_join(icd9_age_diagnosis) -> icd9_diagnosis

icd9_diagnosis %>% 
   filter(!is.na(code_icd9)) -> icd9

feather_file <- file.path(export_path,"ukb_df","icd9.feather")
arrow::write_feather(icd9, feather_file)

# OPCS4 ----

opcs4_codes <- ukbb_df %>% 
   dplyr::select(eid, matches("41272") ) %>% 
   gather("opcs4_diagnosis" , "code_opcs4", - eid) %>%  
   mutate(new = str_remove(opcs4_diagnosis,
                           "operative_procedures_opcs4_f41272_")) 

opcs4_age_diagnosis <- ukbb_df %>% 
   dplyr::select(eid, matches("41282") ) %>% 
   gather("date_version" , "date", - eid) %>%
   mutate(new = str_remove(date_version,
                           "date_of_first_operative_procedure_opcs4_f41282_")) %>% 
   left_join(birth_ukbb_df) %>% 
   mutate(age_opcs4 = time_length(difftime(date, date_of_birth), "years"))

opcs4_codes  %>% 
   left_join(opcs4_age_diagnosis) -> opcs4_diagnosis

opcs4_diagnosis %>% 
   filter(!is.na(code_opcs4)) -> opcs4

feather_file <- file.path(export_path,"ukb_df","opcs4.feather")
arrow::write_feather(opcs4, feather_file)

feather_file <- file.path(export_path,"ukb_df","opcs4.feather")
opcs4_2 <- arrow::read_feather(feather_file) %>% as.data.frame()