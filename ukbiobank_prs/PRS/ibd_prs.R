# 1.0 Load libraries ----

installed.packages("pacman") 
library("pacman")

pacman::p_load(tidyverse,data.table, tictoc, furrr, future,ggpubr, plotly,
               caret, # partition data
               rcompanion, # calculating R2
               pROC,  # calculating and plot AUC 
               plotROC)  # calculating and plot AUC as ggplot   

# 2.0 Importing Files ----

dir_ibd_prs_scores <- "/home/opc/storage/uk_biobank/phenotype_data/prs/ibd"
no_cores <- availableCores() - 1

plan(multisession, workers = no_cores)

ibd_prs_scores <- list.files(dir_ibd_prs_scores , full.names = T) %>%
   future_map(fread) %>%
   reduce(full_join, by= c("FID" = "FID")) %>%
   select(matches("FID"), matches("Pt_1"))

# 3.0 Calculating score sum ----

ibd_prs <-  ibd_prs_scores %>%
   mutate( prs_ibd = rowSums(ibd_prs_scores %>% dplyr::select(- FID), na.rm = TRUE)) %>%
   dplyr::select( FID,prs_ibd) %>%
   rename( eid = `FID`)

prs_pheno <- function(data1, data2){
   
   cases_ids <- data1 %>%
      filter( prev_inc == "prevalent") %>%
      dplyr::select(eid,age,sex_f31_0_0,age_when_attended_assessment_centre_f21003_0_0,
                    genetic_principal_components_f22009_0_1,
                    genetic_principal_components_f22009_0_2,
                    genetic_principal_components_f22009_0_3,
                    genetic_principal_components_f22009_0_4)  %>%
      distinct(.keep_all = TRUE) %>%
      mutate(pheno = 1)
   
   control_ids <- 
      ukbb_df %>% dplyr::select(eid,sex_f31_0_0,ethnic_background_f21000_0_0,
                                outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
                                sex_chromosome_aneuploidy_f22019_0_0,
                                genetic_sex_f22001_0_0,
                                age_when_attended_assessment_centre_f21003_0_0, 
                                genetic_principal_components_f22009_0_1,
                                genetic_principal_components_f22009_0_2,
                                genetic_principal_components_f22009_0_3,
                                genetic_principal_components_f22009_0_4) %>%
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
      dplyr::select(-sex_check) %>%
      anti_join(cases_ids) %>% 
      dplyr::select(eid,age_when_attended_assessment_centre_f21003_0_0,sex_f31_0_0,
                    genetic_principal_components_f22009_0_1,
                    genetic_principal_components_f22009_0_2,
                    genetic_principal_components_f22009_0_3,
                    genetic_principal_components_f22009_0_4)  %>%
      mutate(age = age_when_attended_assessment_centre_f21003_0_0) %>% 
      distinct(.keep_all = TRUE) %>%
      mutate(pheno = 0)
   
   all_ids <- cases_ids %>%
      bind_rows(control_ids)
   
   
   #ukbb_df %>% 
   #   dplyr::select(eid, matches("illnesses_of_mother")) %>%
   #   gather("relationship" , "illness", - eid) %>%
   #   filter(str_detect(illness, 'Breast')) %>% 
   #   dplyr::select(-relationship) %>% 
   #   distinct(eid, .keep_all = TRUE)-> fam_hist
   
   scores_sum_pheno <- data2 %>%
      left_join(all_ids, by = c("eid"="eid")) %>%
      filter(!is.na(pheno)) #%>% 
   #left_join(fam_hist, by = c("id"="eid")) %>% 
   #mutate(ill = ifelse(illness == "Breast cancer", 1, 0 )) %>% 
   #mutate(ill = replace_na(ill,0)) %>% 
   #left_join(ukbb_pc_df %>% dplyr::select(eid, genetic_principal_components_f22009_0_1,
   #                                       genetic_principal_components_f22009_0_2  ),
   #          by = c("id" ="eid"))
   return(scores_sum_pheno)
}

prev_ibd_scores <- prs_pheno(ibd_final,ibd_prs)

predict_auc <- function(data){
   set.seed(123)
   inTrain<- createDataPartition(y=data$pheno, p=0.75, list=FALSE)
   training<-data[inTrain,]
   
   testing<-data[-inTrain,]
   
   glm(pheno ~  
          prs_ibd +
          age_when_attended_assessment_centre_f21003_0_0 +
          sex_f31_0_0 +
          genetic_principal_components_f22009_0_1 +
          genetic_principal_components_f22009_0_2 +
          genetic_principal_components_f22009_0_3 +
          genetic_principal_components_f22009_0_4 , 
       training , family = binomial("logit") , maxit = 100) -> model
   
   predict(model, testing) -> prob
   testing$prob=prob
   # calculating AUC
   g <- roc(pheno ~ prob, data = testing, plot=TRUE,
            print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)  
}

predict_auc(prev_ibd_scores)

genepy_exome_scores %>%  
   inner_join(conversion_table %>% dplyr::select(V1,V7), by = c('IID' = 'V1')) %>% 
   inner_join(prev_ibd_scores %>% select(eid, pheno, prs_ibd, sex_f31_0_0,
                                         age_when_attended_assessment_centre_f21003_0_0,
                                         genetic_principal_components_f22009_0_1,
                                         genetic_principal_components_f22009_0_2,
                                         genetic_principal_components_f22009_0_3,
                                         genetic_principal_components_f22009_0_4)
              , by = c('V7' = 'eid')) %>% 
   mutate(pheno = ifelse(pheno == 1, "yes", "no")) -> prs_gene_score_ibd




laso_select(prs_gene_score_ibd) -> laso_model_ibd

varImp(laso_model_ibd, scale = TRUE)$importance %>% 
   as.data.frame() %>% 
   rownames_to_column() %>%
   arrange(desc(Overall)) %>% 
   filter(Overall != 0) %>% 
   select(1) -> varbs2
varbs2$rowname -> varbs2

