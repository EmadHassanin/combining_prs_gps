# 1.0 Load libraries ----

installed.packages("pacman") 
library("pacman")

pacman::p_load(tidyverse, data.table, tictoc, furrr, future,ggpubr, plotly)

# 2.0 Importing Files ----

dir_scores <- "/home/opc/storage/output/score_313"

no_cores <- availableCores() - 1

plan(multisession, workers = no_cores)

scores <- list.files(dir_scores,pattern = "sscore$", full.names = T) %>%
   future_map(fread) %>%
   reduce(full_join, by= c("#FID" = "#FID")) %>%
   select(matches("#FID"), matches("SCORE"))

# 3.0 calculating score mean ----

scores_mean <-  scores %>%
   mutate( mean =  rowMeans(scores %>% select(- `#FID`), na.rm = TRUE)) %>%
   select(`#FID`,mean) %>%
   rename( id = `#FID`)

cases_ids <- breast_fem_white %>%
   select(eid) %>%
   unique() %>%
   mutate(pheno = 1)

cases_ids <- breast_fem_white %>%
   select(eid) %>%
   unique() %>%
   mutate(pheno = 1)

control_ids <- 
   ukbb_df %>% select(eid,sex_f31_0_0,ethnic_background_f21000_0_0,
                      outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
                      sex_chromosome_aneuploidy_f22019_0_0,
                      genetic_sex_f22001_0_0 ) %>%
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
   select(-sex_check) %>%
   anti_join(cases_ids) %>%
   select(eid)  %>%
   unique() %>%
   mutate(pheno = 0)

all_ids <- cases_ids %>%
   bind_rows(control_ids)




scores_mean_pheno <- scores_mean %>%
   left_join(all_ids, by = c("id"="eid")) %>%
   filter(!is.na(pheno))

# 4.0 exploratoey plots ----

hist(scores_mean_pheno$mean)

scores_mean_pheno %>%
   ggplot(aes(x=as.factor(pheno), y=mean, fill=as.factor(pheno))) +
   geom_boxplot() +
   theme_pubr()

t.test( x= scores_mean_pheno %>% filter(pheno == 1) %>% select(mean),
        y= scores_mean_pheno %>% filter(pheno == 0)%>% select(mean))

scores_mean %>%
   left_join(ukbb_pc_df %>% select(eid, genetic_principal_components_f22009_0_1,
                                   genetic_principal_components_f22009_0_2,
                                   genetic_ethnic_grouping_f22006_0_0),
             by = c("id" ="eid")) %>%
   mutate(eth = case_when(
      genetic_ethnic_grouping_f22006_0_0 == "Caucasian" ~ 1,
      TRUE ~ 0
   )) %>%
   select(mean,eth) -> eth_score

eth_score %>%
   ggplot(aes(x=as.factor(eth), y=mean, fill=as.factor(eth))) +
   geom_boxplot() +
   # geom_point(aes(colour = pheno)) +
   theme_pubr()

t.test( x= eth_score %>% filter(eth == 1) %>% select(mean),
        y= eth_score %>% filter(eth == 0)%>% select(mean))

ukbb_pc_df %>%
   ggplot(aes(genetic_principal_components_f22009_0_1, genetic_principal_components_f22009_0_2,
              color= ethnic_background_f21000_0_0)) +
   geom_point()+
   theme_pubr()

scores_mean_pheno %>%
   left_join(ukbb_pc_df %>% select(eid, genetic_principal_components_f22009_0_1,
                                   genetic_principal_components_f22009_0_2,
                                   age_when_attended_assessment_centre_f21003_0_0  ),
             by = c("id" ="eid")) -> score_pc

# 5.0 regression ----

# association model
scores_mean_pheno %>%
   left_join(ukbb_pc_df %>% select(eid, genetic_principal_components_f22009_0_1,
                                   genetic_principal_components_f22009_0_2,
                                   age_when_attended_assessment_centre_f21003_0_0  ),
             by = c("id" ="eid")) -> score_pc

model <- glm(pheno ~ age_when_attended_assessment_centre_f21003_0_0 + mean + 
                genetic_principal_components_f22009_0_1 +
                genetic_principal_components_f22009_0_2, 
             score_pc , family = "binomial") 

summary(model)
fitted.values(model) %>%  length




# training model

library(caret)

# create test/train data sets

inTrain<- createDataPartition(y=score_pc$pheno, p=0.75, list=FALSE)
training<-score_pc[inTrain,]
testing<-score_pc[-inTrain,]

# fit classification as a model
modelFit<-train(as.factor(pheno) ~ mean + age_when_attended_assessment_centre_f21003_0_0 +
                   genetic_principal_components_f22009_0_1 +
                   genetic_principal_components_f22009_0_2, 
                data=training, method='glm')

predictions<- predict(modelFit, newdata=testing)

glm(pheno ~ mean + age_when_attended_assessment_centre_f21003_0_0 +
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2,
    data=training, family = "binomial") -> model

predict(model, testing) -> prob

testing$prob=prob

# calculating AUC

library(pROC)
g <- roc(pheno ~ prob, data = testing, plot=TRUE,
         print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

library(plotROC)
rocplot <- ggplot(score_pc_tail, aes(m = prob, d = pheno ))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot +
   style_roc(guide = FALSE) +
   geom_rocci()+
   annotate("text", x = .75, y = .25,
            label = paste("AUC =", round(calc_auc(rocplot)$AUC, 2))) +
   theme_pubr()


glm(pheno ~ mean + age_when_attended_assessment_centre_f21003_0_0 +
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2,
    score_pc %>% head(n=300000) %>% filter(!is.na(pheno)),
    family = "binomial") -> model

predict(model, data.frame(score_pc %>% tail(n=100000) %>%
                             filter(!is.na(pheno)))) -> prob

score_pc %>% tail(n=100000) %>%
   filter(!is.na(pheno)) -> score_pc_tail
score_pc_tail$prob=prob

# calculating AUC

library(pROC)
g <- roc(pheno ~ prob, data = score_pc_tail, plot=TRUE,
         print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

library(plotROC)
rocplot <- ggplot(score_pc_tail, aes(m = prob, d = pheno ))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot +
   style_roc(guide = FALSE) +
   geom_rocci()+
   annotate("text", x = .75, y = .25,
            label = paste("AUC =", round(calc_auc(rocplot)$AUC, 2))) +
   theme_pubr()

