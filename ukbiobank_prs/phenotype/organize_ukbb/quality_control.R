pheno_icd_1 %>% 
   mutate( min_age_at_cancer_diagnosis =  apply(pheno_icd_1 %>% select( contains("age_at_cancer_diagnosis")),
                                                1,min, na.rm = TRUE),
           min_age_at_cancer_diagnosis = ifelse(is.infinite(min_age_at_cancer_diagnosis ), 
                                                NA, min_age_at_cancer_diagnosis )) %>% 
   filter(
      genetic_ethnic_grouping_f22006_0_0 == "Caucasian" &
      ethnic_background_f21000_0_0 == "British") %>% #&
      #outliers_for_heterozygosity_or_missing_rate_f22027_0_0 != "Yes" &
      #sex_chromosome_aneuploidy_f22019_0_0 != "Yes"  
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
   distinct(eid) %>% nrow

ukbb_df %>% select( contains("age_at_cancer_diagnosis")) %>% head(n=500) %>% 

 mutate( test =  apply(.,1,min, na.rm = TRUE),
         test = ifelse(is.infinite(test), NA, test)) %>% view

inner_join( ukbb_df  %>% 
               select( col_vars$col_name,  matches("age_at_cancer_diagnosis") ,
                       matches("interpolated_age_of_participant_when_cancer_first")),
            icdcodes_all %>% 
               filter(str_detect(icd9_code , identify_icd$icd9)   | 
                         str_detect(icd10_code ,identify_icd$icd10)),
            by = c( "eid" = "sample" )) %>% 
   
pheno_icd_1_4 %>%
   
   mutate( min_age_at_cancer_diagnosis =  apply(pheno_icd_1_4%>% select( contains("age_at_cancer_diagnosis"),
                                             matches("interpolated_age_of_participant_when_cancer_first")),
                                                1,min, na.rm = TRUE),
           min_age_at_cancer_diagnosis = ifelse(is.infinite(min_age_at_cancer_diagnosis ), 
                                                NA, min_age_at_cancer_diagnosis )) %>% 
   mutate( prev_inc = case_when(
      age_when_attended_assessment_centre_f21003_0_0 > min_age_at_cancer_diagnosis  ~ "prevalent",
      age_when_attended_assessment_centre_f21003_0_0 < min_age_at_cancer_diagnosis  ~ "incident",
      age_when_attended_assessment_centre_f21003_0_0 == min_age_at_cancer_diagnosis ~ "unknown"
   )) %>% 
   filter(sex_f31_0_0 == "Female") %>% 
   filter(
#      genetic_ethnic_grouping_f22006_0_0 == "Caucasian" &
         ethnic_background_f21000_0_0 == "British") %>% 

   inner_join(fam_df %>% select(V1), by = c( "eid" = "V1" )) %>% 
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
  # count(prev_inc) %>% 
    filter(prev_inc == "prevalent") %>% 
   distinct(eid) %>% nrow


pheno_icd_all %>%
   
   mutate( min_age_at_cancer_diagnosis =  apply(pheno_icd_1_4%>% select( contains("age_at_cancer_diagnosis"),
                                                                         matches("interpolated_age_of_participant_when_cancer_first")),
                                                1,min, na.rm = TRUE),
           min_age_at_cancer_diagnosis = ifelse(is.infinite(min_age_at_cancer_diagnosis ), 
                                                NA, min_age_at_cancer_diagnosis )) %>% 
   mutate( prev_inc = case_when(
      age_when_attended_assessment_centre_f21003_0_0 > min_age_at_cancer_diagnosis  ~ "prevalent",
      age_when_attended_assessment_centre_f21003_0_0 < min_age_at_cancer_diagnosis  ~ "incident",
      age_when_attended_assessment_centre_f21003_0_0 == min_age_at_cancer_diagnosis ~ "unknown",
   ), 
   prev_inc = ifelse(is.na(prev_inc ), 
                     "prevalent", prev_inc )) %>%
   filter(sex_f31_0_0 == "Female") %>% 
   filter(
            genetic_ethnic_grouping_f22006_0_0 == "Caucasian" &
      ethnic_background_f21000_0_0 == "British") %>% 
   
   inner_join(fam_df %>% select(V1), by = c( "eid" = "V1" )) %>% 
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
   filter(prev_inc == "prevalent") %>% 
  nrow

pheno_icd_1_4 %>%
   mutate(behaviour = case_when(
      str_detect(behaviour_of_cancer_tumour_f40012_0_0, "Malignant") ~ "Malignant",
      str_detect(behaviour_of_cancer_tumour_f40012_0_0, "Carcinoma") ~ "Carcinoma in situ",
      str_detect(behaviour_of_cancer_tumour_f40012_0_0, "Benign") ~ "Benign",
      TRUE ~ "unspecified")) %>% 
   distinct(eid, .keep_all = T) %>% 
   count(behaviour )
   filter(sex_f31_0_0 == "Female") %>%
   mutate( min_age_at_cancer_diagnosis =  apply(pheno_icd_1_4 %>% 
                                                   filter(sex_f31_0_0 == "Female") %>% 
                                                   select( contains("age_at_cancer_diagnosis"),
                                                                         matches("interpolated_age_of_participant_when_cancer_first")),
                                                1,min, na.rm = TRUE),
           min_age_at_cancer_diagnosis = ifelse(is.infinite(min_age_at_cancer_diagnosis ), 
                                                NA, min_age_at_cancer_diagnosis )) %>% 
   mutate( prev_inc = case_when(
      age_when_attended_assessment_centre_f21003_0_0 > min_age_at_cancer_diagnosis  ~ "prevalent",
      age_when_attended_assessment_centre_f21003_0_0 < min_age_at_cancer_diagnosis  ~ "incident",
      age_when_attended_assessment_centre_f21003_0_0 == min_age_at_cancer_diagnosis ~ "unknown",
   )) %>%  
   inner_join(fam_df %>% select(V1), by = c( "eid" = "V1" )) %>% 
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
   filter(sex_f31_0_0 == "Female") %>% 
   filter(
      genetic_ethnic_grouping_f22006_0_0 == "Caucasian" &
         ethnic_background_f21000_0_0 == "British") %>% 
   distinct(eid, .keep_all = TRUE) %>% 
   count(prev_inc)

   
   
   
   
   pheno_icd_ukbtools %>%
      mutate(behaviour = case_when(
         str_detect(behaviour_of_cancer_tumour_f40012_0_0, "Malignant") ~ "Malignant",
         str_detect(behaviour_of_cancer_tumour_f40012_0_0, "Carcinoma") ~ "Carcinoma in situ",
         str_detect(behaviour_of_cancer_tumour_f40012_0_0, "Benign") ~ "Benign",
         TRUE ~ "unspecified")) %>% 
      mutate( min_age_at_cancer_diagnosis =  apply(pheno_icd_ukbtools%>% select( contains("age_at_cancer_diagnosis"),
                                                                            matches("interpolated_age_of_participant_when_cancer_first")),
                                                   1,min, na.rm = TRUE),
              min_age_at_cancer_diagnosis = ifelse(is.infinite(min_age_at_cancer_diagnosis ), 
                                                   NA, min_age_at_cancer_diagnosis )) %>%
      mutate( prev_inc = case_when(
         age_when_attended_assessment_centre_f21003_0_0 > min_age_at_cancer_diagnosis  ~ "prevalent",
         age_when_attended_assessment_centre_f21003_0_0 < min_age_at_cancer_diagnosis  ~ "incident",
         age_when_attended_assessment_centre_f21003_0_0 == min_age_at_cancer_diagnosis ~ "unknown",
      )) %>% 
      distinct(eid, .keep_all = T) %>% 
      count(behaviour,prev_inc) %>% 
      spread(prev_inc,n)
      
      
   pheno_icd_1_4 %>%
      mutate(behaviour = case_when(
         str_detect(behaviour_of_cancer_tumour_f40012_0_0, "Malignant") ~ "Malignant",
         str_detect(behaviour_of_cancer_tumour_f40012_0_0, "Carcinoma") ~ "Carcinoma in situ",
         str_detect(behaviour_of_cancer_tumour_f40012_0_0, "Benign") ~ "Benign",
         TRUE ~ "unspecified")) %>% 
      mutate( min_date_at_cancer_diagnosis =  apply(pheno_icd_1_4%>% select( contains("date_of_cancer_diagnosis")),
                                                   1,min, na.rm = TRUE),
              min_date_at_cancer_diagnosis  = ifelse(is.infinite(min_date_at_cancer_diagnosis  ), 
                                                   NA, min_date_at_cancer_diagnosis  )) %>%
      mutate( min_date_of_attending_assessment =  apply(pheno_icd_1_4%>% select( contains("date_of_attending_assessment")),
                                                    1,min, na.rm = TRUE),
              min_date_of_attending_assessment = ifelse(is.infinite(min_date_of_attending_assessment ), 
                                                   NA, min_date_of_attending_assessment )) %>% 
      mutate( prev_inc = case_when(
         min_date_of_attending_assessment > min_date_at_cancer_diagnosis ~ "prevalent",
         min_date_of_attending_assessment < min_date_at_cancer_diagnosis  ~ "incident",
         min_date_of_attending_assessment == min_date_at_cancer_diagnosis ~ "unknown",
      )) %>%
      filter(
         genetic_ethnic_grouping_f22006_0_0 == "Caucasian" &
            ethnic_background_f21000_0_0 == "British") %>%
      distinct(eid, .keep_all = T) %>% 
      count(behaviour,prev_inc) %>% 
      spread(prev_inc,n)
   

   
   
   