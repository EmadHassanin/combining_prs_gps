depression <- self_reported %>% filter(code_selfreported == 1286) %>% 
   mutate(age_selfreported = if_else(age_selfreported < 0 , mean(age_selfreported), age_selfreported )) %>% 
   full_join(icd10 %>% filter(str_detect(code_icd10 , "^(F32[0-9]|F33[0-9])")),
             by=c("eid", "new")) %>%
   mutate(age_dp = case_when(
      !is.na(code_icd10) ~ age_icd10,
      !is.na(code_selfreported) ~ age_selfreported
   )) %>%  
   mutate(type_dp = case_when(
      !is.na(code_icd10) ~ "icd10",
      !is.na(code_selfreported) ~ "selfreported")) %>%  
   group_by(eid) %>% 
   slice(which.min(age_dp)) %>% 
   ungroup() %>% 
   #mutate(age =age_selfreported)  %>% 
   select(eid,code_icd10,code_selfreported,age_dp,type_dp) %>% 
   mutate(dp_status = 1) 


bp.disorder <- self_reported %>% filter(code_selfreported == 1291) %>% 
   mutate(age_selfreported = if_else(age_selfreported < 0 , mean(age_selfreported), age_selfreported )) %>% 
   full_join(icd10 %>% filter(str_detect(code_icd10 , "^(F30[0-9]|F31[0-9])")),
             by=c("eid", "new")) %>%
   mutate(age_bp = case_when(
      !is.na(code_icd10) ~ age_icd10,
      !is.na(code_selfreported) ~ age_selfreported
   )) %>%  
   mutate(type_bp = case_when(
      !is.na(code_icd10) ~ "icd10",
      !is.na(code_selfreported) ~ "selfreported")) %>%  
   group_by(eid) %>% 
   slice(which.min(age_bp)) %>% 
   ungroup() %>% 
   #mutate(age =age_selfreported)  %>% 
   select(eid,code_icd10,code_selfreported,age_bp,type_bp) %>% 
   mutate(bp_status = 1) 

schizophrenia <- self_reported %>% filter(code_selfreported == 1289) %>% 
   mutate(age_selfreported = if_else(age_selfreported < 0 , mean(age_selfreported), age_selfreported )) %>% 
   full_join(icd10 %>% filter(str_detect(code_icd10 , "^(F20[0-9])")),
             by=c("eid", "new")) %>%
   mutate(age_schizo = case_when(
      !is.na(code_icd10) ~ age_icd10,
      !is.na(code_selfreported) ~ age_selfreported
   )) %>%  
   mutate(type_schizo = case_when(
      !is.na(code_icd10) ~ "icd10",
      !is.na(code_selfreported) ~ "selfreported")) %>%  
   group_by(eid) %>% 
   slice(which.min(age_schizo)) %>% 
   ungroup() %>% 
   #mutate(age =age_selfreported)  %>% 
   select(eid,code_icd10,code_selfreported,age_schizo,type_schizo) %>% 
   mutate(schizo_status = 1) 

mental.disorders <-  depression %>% 
   full_join(bp.disorder, by = c("eid"="eid")) %>% 
   full_join(schizophrenia, by = c("eid"="eid")) %>% 
   full_join(ukbb_df , by = c("eid"="eid")) %>% 
   mutate(dp_status = replace_na(dp_status, 0),
          bp_status = replace_na(bp_status, 0),
          schizo_status = replace_na(schizo_status, 0)) %>%
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
   select(-sex_check) 

mental.disorders %>% 
   ggplot(aes(as.factor(dp_status),igf1_f30770_0_0)) +
   geom_boxplot()
