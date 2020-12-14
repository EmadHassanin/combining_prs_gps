
covid_prs <- fread("/home/opc/storage/uk_biobank/phenotype_data/extracted_pheno/covid_prs.txt") %>% 
   mutate(severity = as.factor(severity))

covid_prs %>% filter(result == 1) %>% 
   mutate(age_stratify = case_when(
      age_when_attended_assessment_centre_f21003_0_0 <= 50 ~ "young",
      age_when_attended_assessment_centre_f21003_0_0 > 50  & 
         age_when_attended_assessment_centre_f21003_0_0 <= 60~ "middle",
      TRUE ~ "old"
   )) -> covid_prs_test

predict_auc_covid <- function(data){
set.seed(12)
training_1 <- createDataPartition(data$severity, p = 0.75, list=FALSE)

trainData_1 <- data[training_1,] 
testData_1 <- data[-training_1,]

noprsformula<- glm(severity ~ 
                      age_when_attended_assessment_centre_f21003_0_0
                                +sex
                        #         + blood_group
                                +genetic_principal_components_f22009_0_1
                                +genetic_principal_components_f22009_0_2
                                +genetic_principal_components_f22009_0_3
                                +genetic_principal_components_f22009_0_4
                                +genetic_principal_components_f22009_0_5,
                   data=trainData_1, family = binomial)

prsformula<- glm(severity ~  
                    t2d+cad+asthma+stroke+bmi+
                              genetic_principal_components_f22009_0_1+
                              genetic_principal_components_f22009_0_2+
                              genetic_principal_components_f22009_0_3+
                              genetic_principal_components_f22009_0_4+
                              genetic_principal_components_f22009_0_5,
                 data=trainData_1, family = binomial)

allformula <- glm(severity ~   t2d+cad+asthma+stroke+bmi+
                              age_when_attended_assessment_centre_f21003_0_0+
                              sex+
                         #     blood_group+
                              genetic_principal_components_f22009_0_1+
                              genetic_principal_components_f22009_0_2+
                              genetic_principal_components_f22009_0_3+
                              genetic_principal_components_f22009_0_4+
                              genetic_principal_components_f22009_0_5,
                 data=trainData_1, family = binomial)

predict(noprsformula, testData_1) -> prob_a
predict(prsformula, testData_1) -> prob_b
predict(allformula, testData_1) -> prob_c


testData_1$prob_a=prob_a
testData_1$prob_b=prob_b
testData_1$prob_c=prob_c

g_1 <- roc(severity ~ prob_a, data = testData_1, 
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

g_2 <- roc(severity ~ prob_b, data = testData_1,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

g_3 <- roc(severity ~ prob_c, data = testData_1, 
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
c(g_1$auc,g_2$auc,g_3$auc)
#roc.test(g_1,g_3)
}

covid_prs_test %>%
   group_by(age_stratify) %>% 
   nest() %>% 
   mutate(auc =map(data,predict_auc_covid)) -> covid_model_final

covid_model_final %>% unnest(auc)
