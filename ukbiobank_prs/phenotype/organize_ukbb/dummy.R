cases_ids <- breast_fem_indian %>%
   dplyr::select(eid) %>%
   unique() %>%
   mutate(pheno = 1)


control_ids <- 
   ukbb_df %>% dplyr::select(eid,sex_f31_0_0,ethnic_background_f21000_0_0,
                      outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
                      sex_chromosome_aneuploidy_f22019_0_0,
                      genetic_sex_f22001_0_0 ) %>%
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
   dplyr::select(-sex_check) %>%
   anti_join(cases_ids) %>%
   dplyr::select(eid)  %>%
   unique() %>%
   mutate(pheno = 0)

all_ids <- cases_ids %>%
   bind_rows(control_ids)




scores_sum_pheno <- scores_sum %>%
   left_join(all_ids, by = c("id"="eid")) %>%
   filter(!is.na(pheno))




# CAD ----


pheno_cad_icd_ukbtools <- 
   inner_join( ukbb_df  %>% 
                  dplyr::select( col_vars$col_name,  matches("age_at_cancer_diagnosis") ,
                                 matches("interpolated_age_of_participant_when_cancer_first"),
                                 matches("date_of_cancer_diagnosis"),
                                 behaviour_of_cancer_tumour_f40012_0_0,
                                 matches("date_of_attending_assessment_centre")),
               icdcodes_all %>% 
                  filter(str_detect(icd9_code , identify_icd$icd9[2])   | 
                            str_detect(icd10_code ,identify_icd$icd10[2])),
               by = c( "eid" = "sample" ))  %>% 
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
   dplyr::select(-sex_check)

pheno_diabetes_icd_ukbtools <- 
   inner_join( ukbb_df  %>% 
                  dplyr::select( col_vars$col_name,  matches("age_at_cancer_diagnosis") ,
                                 matches("interpolated_age_of_participant_when_cancer_first"),
                                 matches("date_of_cancer_diagnosis"),
                                 behaviour_of_cancer_tumour_f40012_0_0,
                                 matches("date_of_attending_assessment_centre")),
               icdcodes_all %>% 
                  filter( 
                     str_detect(icd10_code ,identify_icd$icd10[3])),
               by = c( "eid" = "sample" ))  %>% 
   #  inner_join(fam_df %>% select(V1), by = c( "eid" = "V1" )) %>% 
   #   full_join(self_reported_non_cancer_diagnosis(ukbb_df, 1248))
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
   dplyr::select(-sex_check)

inner_join(pheno_diabetes_icd_ukbtools,fam_df, by = c( "eid" = "V1" )) 






identify_icd<- 
   tribble(
      ~phenotype  , ~icd10, ~pheno, ~geno,
      "Cardiomyopathy","^(I42[0-9])","1698","1648",
      "Atrioventricular and left bundle-branch bloc","^(I44[0-9])","5430","5269",
      "long QT syndrome / other spec. conduction disorders","^(I458)","87","83",
      "Cardiac arrest","^(I46[0-9])","1428","1378",
      "Ventricular tachycardia","^(I472)","1145","1124",
      "Atrial fibrillation and flutter","^(I48)","19777","19204",
      "Ventricular fibrillation and flutter","^(I490)","455","446",
      "Ventricular premature depolarization","^(I493)","646","633",
      "Heart failure","^(I50[0-9])","7883","7599",
      "Cardiomegaly","^(I517)","3862","3752",
      "Hypertension","^(I10)","98862","95932",
      "Hypertensive heart disease","^(I11[0-9]|I13[0-9])","249","241",
      "Mitral valve disease","^(I34[0-9])","3322","3210",
      "Aortic valve disease","^(I35[0-9])","3643","3531",
      "Infartion/ischemic heart disease","^(I2[1-5])","32115","31185"
   )


inner_join( ukbb_df  %>% 
               dplyr::select( col_vars$col_name),
            icdcodes_all %>% 
               filter( str_detect(icd10_code ,identify_icd$icd10[16])),
            by = c( "eid" = "sample" ))  %>%
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
   dplyr::select(-sex_check) %>%  distinct(eid) %>% nrow() 

fwrite(identify_icd, file = file.path(export_path,
                                      paste('pheno_cardio.txt', sep ="")), 
       sep=",", row.names = FALSE)   