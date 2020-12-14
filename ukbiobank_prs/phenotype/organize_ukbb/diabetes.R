diabetes_all <- self_reported %>% filter(code_selfreported == 1223) %>% 
   mutate(age_selfreported = if_else(age_selfreported < 0 , mean(age_selfreported), age_selfreported )) %>% 
   full_join(icd10 %>% filter(str_detect(code_icd10 , "^(E11[0-9])")),
             by=c("eid", "new")) %>%
   mutate(age = case_when(
      !is.na(code_icd10) ~ age_icd10,
      !is.na(code_selfreported) ~ age_selfreported
   )) %>%  
   mutate(type = case_when(
      !is.na(code_icd10) ~ "icd10",
      !is.na(code_selfreported) ~ "selfreported")) %>%  
   group_by(eid) %>% 
   slice(which.min(age)) %>% 
   ungroup() %>% 
   #mutate(age =age_selfreported)  %>% 
   select(eid,code_icd10,code_selfreported,age,type) %>% 
   mutate(pheno = 1) 

diabetes_final <- 
   ukbb_df %>% select(
      eid,sex_f31_0_0, ethnic_background_f21000_0_0, 
      matches("f22009"), sex_chromosome_aneuploidy_f22019_0_0,
      outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
      genetic_sex_f22001_0_0,
      age_when_attended_assessment_centre_f21003_0_0) %>% 
   inner_join(diabetes_all) %>% 
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
      age_when_attended_assessment_centre_f21003_0_0 >= age  ~ "prevalent",
      age_when_attended_assessment_centre_f21003_0_0 < age  ~ "incident",
      TRUE ~ "unknown"
   )) %>% 
   select(
      eid,  age, ethnic_background_f21000_0_0, type,age_when_attended_assessment_centre_f21003_0_0,
      pheno,prev_inc,sex_f31_0_0, matches("f22009")) 

