af_all <- self_reported %>% filter(code_selfreported == 1471 | 
                                      code_selfreported == 1483) %>% 
   mutate(age_selfreported = if_else(age_selfreported < 0 , mean(age_selfreported), age_selfreported )) %>% 
   full_join(icd9 %>% filter(str_detect(code_icd9 , "^(4273)")),
             by=c("eid", "new")) %>% 
   full_join(icd10 %>% filter(str_detect(code_icd10 , "^(I48[0-9])")),
             by=c("eid", "new")) %>%
   full_join(opcs4 %>% filter(str_detect(code_opcs4 , 
                                         "^(K571|K62[1-4])")),
             by=c("eid", "new")) %>%
   mutate(age = case_when(
      !is.na(code_icd10) ~ age_icd10,
      !is.na(code_opcs4) ~ age_opcs4,
      !is.na(code_selfreported) ~ age_selfreported,
      !is.na(code_icd9) ~ age_icd9
   )) %>%  
   mutate(type = case_when(
      !is.na(code_icd10) ~ "icd10",
      !is.na(code_opcs4) ~ "opcs4",
      !is.na(code_selfreported) ~ "selfreported",
      !is.na(code_icd9) ~ "icd9")) %>%  
   group_by(eid) %>% 
   slice(which.min(age)) %>% 
   ungroup() %>% 
   #mutate(age =age_selfreported)  %>% 
   select(eid,code_icd10,code_icd9,code_opcs4,code_selfreported,age,type) %>% 
   mutate(pheno = 1) 

af_final <- 
   ukbb_df %>% select(
      eid,sex_f31_0_0, ethnic_background_f21000_0_0, 
      matches("f22009"), sex_chromosome_aneuploidy_f22019_0_0,
      outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
      genetic_sex_f22001_0_0,
      age_when_attended_assessment_centre_f21003_0_0) %>% 
   inner_join(af_all) %>% 
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

