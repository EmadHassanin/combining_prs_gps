# 1.0 Load libraries ----
pacman::p_load(tidyverse     # tidy daya
               , ukbtools    # workhorse package for organizing ukbb
               , tictoc      # checking time 
               , furrr       # functional programing with map()
               , future      # parallel computing
               , arrow       # apache arrow provides powerful package for sharing data
               , data.table  # exporting data as .txt
               , ggpubr      # add themes and visualization prop. to ggplot
               ,plotly       # interactive plots
               ,caret        # partition data
               #,rcompanion   # calculating R2
               ,pROC         # calculating and plot AUC 
               ,plotROC      # calculating and plot ggplot
               ,gridExtra    # add table to ggplot
)    

# 2.0 functions

self_reported_cancer_diagnosis <- function(data, coding){
data %>% 
   select(eid, matches("cancer_code_selfreported") ) %>% 
   gather("cancer_code_selfreported_version" , "code", - eid) %>% 
   select(-cancer_code_selfreported_version)   %>% 
   filter(code == coding) %>% 
   distinct(eid, code)  
}

register_cancer_diagnosis <- function(data, coding){
   data %>% 
      select(eid, matches("type_of_cancer") ) %>% 
      gather("type_of_cancer_version" , "icd_code", - eid) %>%
      select(-type_of_cancer_version)   %>% 
      filter(str_detect(icd_code , coding)) %>%
 #     filter(!is.na(icd_code)) %>% 
      distinct(eid, icd_code)  
}



self_reported_non_cancer_diagnosis <- function(data, coding){
   data %>% 
      dplyr::select(eid, matches("noncancer_illness_code_selfreported") ) %>% 
      gather("non_cancer_code_selfreported_version" , "code", - eid) %>% 
      dplyr::select(-non_cancer_code_selfreported_version)   %>% 
      filter(code == coding) %>% 
      distinct(eid, code)  
}



predict_auc <- function(data){
set.seed(12)
inTrain<- createDataPartition(y=data$pheno, p=0.75, list=FALSE)
training<-data[inTrain,]

testing<-data[-inTrain,]

glm(pheno ~  
       sum + 
       age +
  #    ill +
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2, 
    training , family = "binomial") -> model

predict(model, testing) -> prob
testing$prob=prob
# calculating AUC
g <- roc(pheno ~ prob, data = testing, plot=TRUE,
         print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)  
}



#tribble(
#   ~option, ~`313_snps`,~khera_snps_prsice,
#   "all", 0.63 , 0.63,
#   "incindet", 0.684, 0.69,
#   "prevalent", 0.712, 0.711,
#   "prevalent + genepy", 0.732, 0.725,
#   "prevalent + incident",0.828, 0.828
#) 

prs_pheno <- function(data, type){

cases_ids <- breast_fem_white  %>%
   #filter( prev_inc == "incident") %>%
   filter( prev_inc == type) %>%
   dplyr::select(eid,min_age_at_cancer_diagnosis)  %>%
   rename(age = min_age_at_cancer_diagnosis) %>% 
   distinct(.keep_all = TRUE) %>%
   mutate(pheno = 1)

control_ids <- 
   ukbb_df %>% dplyr::select(eid,sex_f31_0_0,ethnic_background_f21000_0_0,
                             outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
                             sex_chromosome_aneuploidy_f22019_0_0,
                             genetic_sex_f22001_0_0,
                             age_when_attended_assessment_centre_f21003_0_0) %>%
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
   dplyr::select(-sex_check) %>%
   anti_join(cases_ids) %>% 
   dplyr::select(eid,age_when_attended_assessment_centre_f21003_0_0)  %>%
   rename(age = age_when_attended_assessment_centre_f21003_0_0) %>% 
   distinct(.keep_all = TRUE) %>%
   mutate(pheno = 0)

all_ids <- cases_ids %>%
   bind_rows(control_ids)





ukbb_df %>% 
   dplyr::select(eid, matches("illnesses_of_mother")) %>%
   gather("relationship" , "illness", - eid) %>%
   filter(str_detect(illness, 'Breast')) %>% 
   dplyr::select(-relationship) %>% 
   distinct(eid, .keep_all = TRUE)-> fam_hist

scores_sum_pheno <- data %>%
   left_join(all_ids, by = c("id"="eid")) %>%
   filter(!is.na(pheno)) %>% 
   left_join(fam_hist, by = c("id"="eid")) %>% 
   mutate(ill = ifelse(illness == "Breast cancer", 1, 0 )) %>% 
   mutate(ill = replace_na(ill,0)) %>% 
   left_join(ukbb_pc_df %>% dplyr::select(eid, genetic_principal_components_f22009_0_1,
                                          genetic_principal_components_f22009_0_2  ),
             by = c("id" ="eid"))
return(scores_sum_pheno)
}

prs_pheno_all <- function(data){
   
   cases_ids <- breast_fem_white  %>%
      dplyr::select(eid,min_age_at_cancer_diagnosis)  %>%
      rename(age = min_age_at_cancer_diagnosis) %>% 
      distinct(.keep_all = TRUE) %>%
      mutate(pheno = 1)
   
   control_ids <- 
      ukbb_df %>% dplyr::select(eid,sex_f31_0_0,ethnic_background_f21000_0_0,
                                outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
                                sex_chromosome_aneuploidy_f22019_0_0,
                                genetic_sex_f22001_0_0,
                                age_when_attended_assessment_centre_f21003_0_0) %>%
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
      dplyr::select(-sex_check) %>%
      anti_join(cases_ids) %>% 
      dplyr::select(eid,age_when_attended_assessment_centre_f21003_0_0)  %>%
      rename(age = age_when_attended_assessment_centre_f21003_0_0) %>% 
      distinct(.keep_all = TRUE) %>%
      mutate(pheno = 0)
   
   all_ids <- cases_ids %>%
      bind_rows(control_ids)
   
   
   
   
   
   ukbb_df %>% 
      dplyr::select(eid, matches("illnesses_of_mother")) %>%
      gather("relationship" , "illness", - eid) %>%
      filter(str_detect(illness, 'Breast')) %>% 
      dplyr::select(-relationship) %>% 
      distinct(eid, .keep_all = TRUE)-> fam_hist
   
   scores_sum_pheno <- data %>%
      left_join(all_ids, by = c("id"="eid")) %>%
      filter(!is.na(pheno)) %>% 
      left_join(fam_hist, by = c("id"="eid")) %>% 
      mutate(ill = ifelse(illness == "Breast cancer", 1, 0 )) %>% 
      mutate(ill = replace_na(ill,0)) %>% 
      left_join(ukbb_pc_df %>% dplyr::select(eid, genetic_principal_components_f22009_0_1,
                                             genetic_principal_components_f22009_0_2  ),
                by = c("id" ="eid"))
   return(scores_sum_pheno)
}

prs_pheno_indian <- function(data,type){
   
   cases_ids <- breast_fem_indian %>%
      #filter( prev_inc == "incident") %>%
      filter( prev_inc == type) %>%
      dplyr::select(eid,min_age_at_cancer_diagnosis)  %>%
      rename(age = min_age_at_cancer_diagnosis) %>% 
      distinct(.keep_all = TRUE) %>%
      mutate(pheno = 1)
   
   control_ids <- 
      ukbb_df %>% dplyr::select(eid,sex_f31_0_0,ethnic_background_f21000_0_0,
                                outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
                                sex_chromosome_aneuploidy_f22019_0_0,
                                genetic_sex_f22001_0_0,
                                age_when_attended_assessment_centre_f21003_0_0) %>%
      filter(sex_f31_0_0 == "Female" & ethnic_background_f21000_0_0 =="Indian" ) %>%
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
      dplyr::select(eid,age_when_attended_assessment_centre_f21003_0_0)  %>%
      rename(age = age_when_attended_assessment_centre_f21003_0_0) %>% 
      distinct(.keep_all = TRUE) %>%
      mutate(pheno = 0)
   
   all_ids <- cases_ids %>%
      bind_rows(control_ids)
   
   
   
   
   
   ukbb_df %>% 
      dplyr::select(eid, matches("illnesses_of_mother")) %>%
      gather("relationship" , "illness", - eid) %>%
      filter(str_detect(illness, 'Breast')) %>% 
      dplyr::select(-relationship) %>% 
      distinct(eid, .keep_all = TRUE)-> fam_hist
   
   scores_sum_pheno <- data %>%
      left_join(all_ids, by = c("id"="eid")) %>%
      filter(!is.na(pheno)) %>% 
      left_join(fam_hist, by = c("id"="eid")) %>% 
      mutate(ill = ifelse(illness == "Breast cancer", 1, 0 )) %>% 
      mutate(ill = replace_na(ill,0)) %>% 
      left_join(ukbb_pc_df %>% dplyr::select(eid, genetic_principal_components_f22009_0_1,
                                             genetic_principal_components_f22009_0_2  ),
                by = c("id" ="eid"))
   return(scores_sum_pheno)
}

laso_select <- function(data){

set.seed(1)
training <- createDataPartition(data$pheno, p = 0.75, list=FALSE)

trainData <- data[training,] 
testData <- data[-training,]


fitControl <- trainControl(method = "cv",
                           number = 5,
#                           repeats = 5,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary
)

set.seed(1)
myGrid <- expand.grid(
   alpha = 1,
   lambda = 10^seq(1, -7, length = 10))

library(doParallel)
cl <- makePSOCKcluster(2)
registerDoParallel(cl)
lassoModel20 <- train(x=trainData %>% 
                         dplyr::select(-pheno, -IID, -prs_ibd, 
                         -sex_f31_0_0,
                         -age_when_attended_assessment_centre_f21003_0_0,
                         -genetic_principal_components_f22009_0_1,
                         -genetic_principal_components_f22009_0_2,
                         -genetic_principal_components_f22009_0_3,
                         -genetic_principal_components_f22009_0_4), y=trainData$pheno,
                     method = "glmnet",  metric="ROC",
                     trControl = fitControl, verbose=FALSE, tuneGrid = myGrid ,
                     preProcess = c("zv","medianImpute","center","scale") )
}


## plot gps+prs for Breast cancer

integrate_gps_prs <- function(training,testing){
   glm(as.factor(pheno)~  PRS + exome_score_sum + age
    + genetic_principal_components_f22009_0_1 
    + genetic_principal_components_f22009_0_2,
    
    data=training, family = binomial) -> model_1



glm(as.factor(pheno) ~ PRS +  age +
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2,
    data=training, family = binomial) -> model_2

glm(as.factor(pheno) ~ exome_score_sum  + age +
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2,
    data=training, family = "binomial") -> model_3



predict(model_1, testing) -> prob_1
predict(model_2, testing) -> prob_2
predict(model_3, testing) -> prob_3


testing$prob_1=prob_1
testing$prob_2=prob_2
testing$prob_3=prob_3


longtest <- melt_roc(testing, "pheno", c("prob_1", "prob_2","prob_3")) %>% 
   mutate(name = case_when(
      name == "prob_1" ~ "gps + prs", 
      name == "prob_2" ~"prs",
      name == "prob_3" ~"gps"))

model <- unique(longtest$name)
model_info <- data.frame(model,
                         group = rank(model))

ggplot(longtest, aes(d = D, m = M, color = name)) + 
   geom_roc(n.cuts=20,labels=FALSE) ->p

left_join(model_info, calc_auc(p)) %>%
   select(-group, -PANEL) %>%
   arrange(desc(AUC)) %>% 
   mutate( AUC = round(AUC,3)) -> aucs



p +
   annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()), 
                     xmin=0.5, xmax=1, ymin=0, ymax=0.5) + 
   theme_pubr()

}

gps_prs <- function(training,testing){
   glm(as.factor(pheno) ~ PRS +
          genetic_principal_components_f22009_0_1 +
          genetic_principal_components_f22009_0_2,
       data=training, family = "binomial") -> model_4
   
   glm(as.factor(pheno) ~ exome_score_sum +
          genetic_principal_components_f22009_0_1 +
          genetic_principal_components_f22009_0_2,
       data=training, family = "binomial") -> model_5
 
   predict(model_4, testing) -> prob_4
   predict(model_5, testing) -> prob_5  
   
   testing$prob_4=prob_4
   testing$prob_5=prob_5
   longtest2 <- melt_roc(testing, "pheno", c("prob_4", "prob_5")) %>% 
      mutate(name = case_when(
         name == "prob_4" ~ "prs only", 
         name == "prob_5" ~"gps only"))
   
   model <- unique(longtest2$name)
   model_info <- data.frame(model,
                            group = rank(model))
   
   ggplot(longtest2, aes(d = D, m = M, color = name)) + 
      geom_roc(n.cuts=20,labels=FALSE) ->p2
   
   left_join(model_info, calc_auc(p2)) %>%
      select(-group, -PANEL) %>%
      arrange(desc(AUC)) %>% 
      mutate( AUC = round(AUC,3)) -> aucs  
   p2 +
      annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()), 
                        xmin=0.5, xmax=1, ymin=0, ymax=0.5) + 
      theme_pubr()
}

## plot gps+prs for CAD

integrate_cad_gps_prs <- function(training,testing){
glm(pheno ~  prs_cad + exome_score_sum + age_when_attended_assessment_centre_f21003_0_0
    +sex_f31_0_0
    + genetic_principal_components_f22009_0_1 
    + genetic_principal_components_f22009_0_2
    +genetic_principal_components_f22009_0_3
    +genetic_principal_components_f22009_0_4,
    
    data=training, family = binomial) -> model_1



glm(pheno ~ prs_cad + age_when_attended_assessment_centre_f21003_0_0
    +sex_f31_0_0
    + genetic_principal_components_f22009_0_1 
    + genetic_principal_components_f22009_0_2
    +genetic_principal_components_f22009_0_3
    +genetic_principal_components_f22009_0_4,
    data=training, family = binomial) -> model_2

glm(pheno ~  exome_score_sum + age_when_attended_assessment_centre_f21003_0_0
    +sex_f31_0_0
    + genetic_principal_components_f22009_0_1 
    + genetic_principal_components_f22009_0_2
    +genetic_principal_components_f22009_0_3
    +genetic_principal_components_f22009_0_4,
    data=training, family = "binomial") -> model_3



predict(model_1, testing) -> prob_1
predict(model_2, testing) -> prob_2
predict(model_3, testing) -> prob_3


testing$prob_1=prob_1
testing$prob_2=prob_2
testing$prob_3=prob_3


longtest <- melt_roc(testing, "pheno", c("prob_1", "prob_2","prob_3")) %>% 
   mutate(name = case_when(
      name == "prob_1" ~ "gps + prs", 
      name == "prob_2" ~"prs",
      name == "prob_3" ~"gps"))

model <- unique(longtest$name)
model_info <- data.frame(model,
                         group = rank(model))

ggplot(longtest, aes(d = D, m = M, color = name)) + 
   geom_roc(n.cuts=20,labels=FALSE) ->p

left_join(model_info, calc_auc(p)) %>%
   select(-group, -PANEL) %>%
   arrange(desc(AUC)) %>% 
   mutate( AUC = round(AUC,3)) -> aucs



p +
   annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()), 
                     xmin=0.5, xmax=1, ymin=0, ymax=0.5) + 
   theme_pubr()
}


cad_gps_prs <- function(training,testing){
   glm(pheno ~ prs_cad 
       + genetic_principal_components_f22009_0_1 
       + genetic_principal_components_f22009_0_2
       + genetic_principal_components_f22009_0_3
       + genetic_principal_components_f22009_0_4,
       data=training, family = "binomial") -> model_4
   
   glm(pheno ~ exome_score_sum 
       + genetic_principal_components_f22009_0_1 
       + genetic_principal_components_f22009_0_2
       + genetic_principal_components_f22009_0_3
       + genetic_principal_components_f22009_0_4,
       data=training, family = "binomial") -> model_5
   
   predict(model_4, testing) -> prob_4
   predict(model_5, testing) -> prob_5  
   
   testing$prob_4=prob_4
   testing$prob_5=prob_5
   longtest2 <- melt_roc(testing, "pheno", c("prob_4", "prob_5")) %>% 
      mutate(name = case_when(
         name == "prob_4" ~ "prs only", 
         name == "prob_5" ~"gps only"))
   
   model <- unique(longtest2$name)
   model_info <- data.frame(model,
                            group = rank(model))
   
   ggplot(longtest2, aes(d = D, m = M, color = name)) + 
      geom_roc(n.cuts=20,labels=FALSE) ->p2
   
   left_join(model_info, calc_auc(p2)) %>%
      select(-group, -PANEL) %>%
      arrange(desc(AUC)) %>% 
      mutate( AUC = round(AUC,3)) -> aucs  
   p2 +
      annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()), 
                        xmin=0.5, xmax=1, ymin=0, ymax=0.5) + 
      theme_pubr()
}

## plot gps+prs for diabetes

integrate_diabetes_gps_prs <- function(training,testing){
   glm(as.factor(pheno) ~  prs_diabetes + exome_score_sum 
       + age_when_attended_assessment_centre_f21003_0_0
       +sex_f31_0_0
       + genetic_principal_components_f22009_0_1 
       + genetic_principal_components_f22009_0_2
       +genetic_principal_components_f22009_0_3
       +genetic_principal_components_f22009_0_4,
       
       data=training, family = binomial) -> model_1
   
   
   
   glm(as.factor(pheno) ~ prs_diabetes
       + age_when_attended_assessment_centre_f21003_0_0
       +sex_f31_0_0
       + genetic_principal_components_f22009_0_1 
       + genetic_principal_components_f22009_0_2
       +genetic_principal_components_f22009_0_3
       +genetic_principal_components_f22009_0_4,
       data=training, family = binomial) -> model_2
   
   glm(as.factor(pheno) ~  exome_score_sum 
       + age_when_attended_assessment_centre_f21003_0_0
       +sex_f31_0_0
       + genetic_principal_components_f22009_0_1 
       + genetic_principal_components_f22009_0_2
       +genetic_principal_components_f22009_0_3
       +genetic_principal_components_f22009_0_4,
       data=training, family = "binomial") -> model_3
   
   
   
   predict(model_1, testing) -> prob_1
   predict(model_2, testing) -> prob_2
   predict(model_3, testing) -> prob_3
   
   
   testing$prob_1=prob_1
   testing$prob_2=prob_2
   testing$prob_3=prob_3
   
   
   longtest <- melt_roc(testing, "pheno", c("prob_1", "prob_2","prob_3")) %>% 
      mutate(name = case_when(
         name == "prob_1" ~ "gps + prs", 
         name == "prob_2" ~"prs",
         name == "prob_3" ~"gps"))
   
   model <- unique(longtest$name)
   model_info <- data.frame(model,
                            group = rank(model))
   
   ggplot(longtest, aes(d = D, m = M, color = name)) + 
      geom_roc(n.cuts=20,labels=FALSE) ->p
   
   left_join(model_info, calc_auc(p)) %>%
      select(-group, -PANEL) %>%
      arrange(desc(AUC)) %>% 
      mutate( AUC = round(AUC,3)) -> aucs
   
   
   
   p +
      annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()), 
                        xmin=0.5, xmax=1, ymin=0, ymax=0.5) + 
      theme_pubr()
}


diabetes_gps_prs <- function(training,testing){
   glm(as.factor(pheno) ~ prs_diabetes 
       + genetic_principal_components_f22009_0_1 
       + genetic_principal_components_f22009_0_2
       + genetic_principal_components_f22009_0_3
       + genetic_principal_components_f22009_0_4,
       data=training, family = "binomial") -> model_4
   
   glm(as.factor(pheno)~ exome_score_sum 
       + genetic_principal_components_f22009_0_1 
       + genetic_principal_components_f22009_0_2
       + genetic_principal_components_f22009_0_3
       + genetic_principal_components_f22009_0_4,
       data=training, family = "binomial") -> model_5
   
   predict(model_4, testing) -> prob_4
   predict(model_5, testing) -> prob_5  
   
   testing$prob_4=prob_4
   testing$prob_5=prob_5
   longtest2 <- melt_roc(testing, "pheno", c("prob_4", "prob_5")) %>% 
      mutate(name = case_when(
         name == "prob_4" ~ "prs only", 
         name == "prob_5" ~"gps only"))
   
   model <- unique(longtest2$name)
   model_info <- data.frame(model,
                            group = rank(model))
   
   ggplot(longtest2, aes(d = D, m = M, color = name)) + 
      geom_roc(n.cuts=20,labels=FALSE) ->p2
   
   left_join(model_info, calc_auc(p2)) %>%
      select(-group, -PANEL) %>%
      arrange(desc(AUC)) %>% 
      mutate( AUC = round(AUC,3)) -> aucs  
   p2 +
      annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()), 
                        xmin=0.5, xmax=1, ymin=0, ymax=0.5) + 
      theme_pubr()
}
