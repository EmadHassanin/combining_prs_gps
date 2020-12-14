# 1.0 Load libraries ----
pacman::p_load(tidyverse     # tidy daya
               , ukbtools    # workhorse package for organizing ukbb
               , tictoc      # checking time 
               , furrr       # functional programing with map()
               , future      # parallel computing
               , arrow       # apache arrow provides powerful package for sharing data
               , data.table  # exporting data as .txt
               , ggpubr      # add themes and visualization prop. to ggplot
               , plotly       # interactive plots
               , caret        # partition data
               , rcompanion   # calculating R2
               , pROC         # calculating and plot AUC 
               , plotROC      # calculating and plot ggplot
               , gridExtra    # add table to ggplot
)    

# 2.0 import files ----
source(here::here("phenotype","organize_ukbb","change_me.R"))
source(here::here("phenotype","organize_ukbb","test.R"))
source(here::here("phenotype","organize_ukbb","organize_ukbb.R"))
source(here::here("PRS","import_prs_scores.R"))

# 3.0 build model ----

prs_gene_score2 <-  prs_gene_score %>% 
   dplyr::select(-illness,-ill) %>% 
 #  select(IID,pheno,age,sum
 #         ,PALB2,BRCA2,BRCA1,CHEK2
 #         ) %>% 
   mutate(pheno = ifelse(pheno == 1 , "yes" , "no"),
          pheno = as.factor(pheno))

#prs_gene_score4 <- nearZeroVar(prs_gene_score2)
#prs_gene_score5 <- prs_gene_score2[, -prs_gene_score4] %>% mutate(pheno = prs_gene_score2$pheno)

prs_gene_score3 <-  prs_gene_score %>% 
   #select(-illness) %>% 
     select(IID,pheno,age,sum
            ,3:30
            ) %>% 
   mutate(pheno = ifelse(pheno == 1 , "yes" , "no"),
          pheno = as.factor(pheno))

set.seed(1)
training <- createDataPartition(prs_gene_score2$pheno, p = 0.75, list=FALSE)

trainData <- prs_gene_score2[training,] 
testData <- prs_gene_score2[-training,]

set.seed(1)
training2 <- createDataPartition(prs_gene_score3$pheno, p = 0.75, list=FALSE)
trainData2 <- prs_gene_score3[training2,] 
testData2 <- prs_gene_score3[-training2,]

fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary
                           )

# gbm ----
#library(doParallel)
#cl <- makePSOCKcluster(10)
#registerDoParallel(cl)
#gbmModel <- train(x=trainData2 %>% select(-pheno), y=trainData2$pheno, method = "gbm",
#                  metric="ROC", trControl = fitControl, verbose=FALSE, tuneLength=5)

#pred.gbmModel <- as.vector(predict(gbmModel, newdata=testData2, type="prob")[,"yes"])


#roc.gbmModel <- pROC::roc(testData2$pheno, pred.gbmModel)

#auc.gbmModel <- pROC::auc(roc.gbmModel)

#library(doParallel)
#cl <- makePSOCKcluster(10)
#registerDoParallel(cl)
# glmboost
#glmboostModel <- train(x=trainData2 %>% select(-pheno), y=trainData2$pheno, 
#                       method = "glmStepAIC", trControl = fitControl, family = "binomial")

#pred.glmboostModel <- as.vector(predict(glmboostModel, newdata=testData2, type="prob")[,"yes"])


#roc.glmboostModel <- pROC::roc(testData2$pheno, pred.glmboostModel)

#auc.glmboostModel <- pROC::auc(roc.glmboostModel)

# glmnet (lasso) ----
#lambda <- 10^seq(-3, -7, length = 10)
set.seed(1)
myGrid <- expand.grid(
   alpha = 1,
   lambda = lambda <- 10^seq(1, -7, length = 10))


#10^seq(3, -10, length = 10)
library(doParallel)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)
lassoModel <- train(x=trainData %>% dplyr::select(-pheno, -IID), y=trainData$pheno,
                    method = "glmnet",  metric="ROC",
                    trControl = fitControl, verbose=FALSE, tuneGrid = myGrid ,
                    preProcess = c("zv","center","scale") )

lassoModel2 <- train(x=trainData2 %>% dplyr::select(-pheno, -IID), y=trainData$pheno,
                    method = "glmnet",  metric="ROC",
                    trControl = fitControl, verbose=FALSE, tuneGrid = myGrid ,
                    preProcess = c("nzv","zv","center","scale") )

lassoModel3 <- train(x=trainData2 %>% dplyr::select(-pheno, -IID, -sum,-age), y=trainData2$pheno,
                     method = "glmnet",  metric="ROC",
                     trControl = fitControl, verbose=FALSE, tuneGrid = myGrid ,
                     preProcess = c("zv","medianImpute","center","scale") )

lassoModel4 <- train(x=trainData %>% dplyr::select(-pheno, -IID, -sum, -age), y=trainData$pheno,
                     method = "glmnet",  metric="ROC",
                     trControl = fitControl, verbose=FALSE, tuneGrid = myGrid ,
                     preProcess = c("nzv","zv","medianImpute","center","scale") )

lassoModel5 <- train(x=training%>% dplyr::select(-pheno, -IID), y=training$pheno,
                     method = "glmnet",  metric="ROC",
                     trControl = fitControl, verbose=FALSE, tuneGrid = myGrid ,
                     preProcess = c("center","scale") )

pred.lassoModel <- as.vector(predict(lassoModel5, newdata=testing,s=0.0003593814,  type="prob")[,"yes"])

roc.lassoModel <- pROC::roc(testing$pheno, pred.lassoModel)

auc.lassoModel <- pROC::auc(roc.lassoModel)

stopCluster(cl)

varImp(lassoModel3, scale = TRUE)$importance %>% 
        as.data.frame() %>% 
    rownames_to_column() %>%
     arrange(desc(Overall)) %>% 
   filter(Overall != 0) %>% 
      # head(n=20) %>% 
 #  filter(rowname != "ill") %>% 
   select(1) -> varbs2
varbs2$rowname -> varbs2
varbs$rowname -> varbs

prs_gene_score %>% 
   dplyr::select(imp_varbs,pheno,IID,genetic_principal_components_f22009_0_1, 
          genetic_principal_components_f22009_0_2,sum,age)  -> score_interest

#delete_score_interest2 <- nearZeroVar(score_interest)
#score_interest <- score_interest[, -delete_score_interest2] %>% mutate(pheno = score_interest$pheno) %>% 
score_interest <-   score_interest %>% 
   mutate(pheno = ifelse(pheno == 1 , "yes" , "no"),
          pheno = as.factor(pheno)) 


score_interest <- score_interest %>% 
   mutate( exome_score_sum =  rowSums(score_interest %>% dplyr::select(-IID,-age,-sum,-pheno,
                                                                       -genetic_principal_components_f22009_0_1,
                                                   -genetic_principal_components_f22009_0_2), 
                                      na.rm = TRUE))

set.seed(12)
inTrain<- createDataPartition(y=score_interest$pheno, p=0.75, list=FALSE)
training<-score_interest[inTrain,]

testing<-score_interest[-inTrain,]

control_genes <- score_interest %>% 
   filter(pheno == "no") %>% dplyr::select(exome_score_sum)

cases_genes <-score_interest  %>% 
   filter(pheno == "yes") %>% dplyr::select(exome_score_sum)

t.test(control_genes$exome_score_sum,cases_genes$exome_score_sum)

glm(pheno ~  sum + exome_score_sum + age
    + genetic_principal_components_f22009_0_1 
    + genetic_principal_components_f22009_0_2,
    
    data=training, family = binomial) -> model_1

#model_1 %>% stepAIC(direction='both',trace = FALSE) -> model_2

glm(pheno ~ sum +  age +
          genetic_principal_components_f22009_0_1 +
         genetic_principal_components_f22009_0_2,
    data=training, family = binomial) -> model_2

glm(pheno ~ exome_score_sum  + age +
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2,
    data=training, family = "binomial") -> model_3

glm(pheno ~ prs +
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2,
    data=training, family = "binomial") -> model_4

glm(pheno ~ exome_score_sum +
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2,
    data=training, family = "binomial") -> model_5

predict(model_1, testing) -> prob_1
predict(model_2, testing) -> prob_2
predict(model_3, testing) -> prob_3
predict(model_4, testing) -> prob_4
predict(model_5, testing) -> prob_5

testing$prob_1=prob_1
testing$prob_2=prob_2
testing$prob_3=prob_3
testing$prob_4=prob_4
testing$prob_5=prob_5

g_1 <- roc(pheno ~ prob_1, data = testing, plot=TRUE,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

g_2 <- roc(pheno ~ prob_2, data = testing, plot=TRUE,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

g_3 <- roc(pheno ~ prob_3, data = testing, plot=TRUE,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

g_4<- roc(pheno ~ prob_4, data = testing, plot=TRUE,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

g_5 <- roc(pheno ~ prob_5, data = testing, plot=TRUE,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

roc.test(g_1,g_2)

longtest <- melt_roc(testing, "pheno", c("prob_1", "prob_2","prob_3")) %>% 
   mutate(name = case_when(
      name == "prob_1" ~ "gps + prs", 
      name == "prob_2" ~"prs ",
      name == "prob_3" ~"gps ",
      name == "prob_4" ~"prs wo cov",
      name == "prob_5" ~"gps wo cov "))

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

fit <- train(pheno ~  sum + exome_score_sum + age
             + genetic_principal_components_f22009_0_1 
             + genetic_principal_components_f22009_0_2,
             
             data=training, 
             method    = "glm"    ,
             family    = binomial ,
             trControl = fitControl)

pred.lassoModel <- as.vector(predict(fit, newdata=testing , type="prob")[,"yes"])

roc.lassoModel <- pROC::roc(testing$pheno, pred.lassoModel)

auc.lassoModel <- pROC::auc(roc.lassoModel)


laso_select(prs_gene_score2) -> lasomodel_1

varImp(lasomodel_1, scale = TRUE)$importance %>% 
   as.data.frame() %>% 
   rownames_to_column() %>%
   arrange(desc(Overall)) %>% 
   filter(Overall != 0) %>% 
   select(1) ->  imp_varbs
imp_varbs$rowname -> imp_varbs
saveRDS(imp_varbs, file = "/home/opc/storage/output/output/prs_genepy/gps_genes_breast_white.rds")
imp_varbs <- readRDS(file =  "/home/opc/storage/output/output/prs_genepy/gps_genes_breast_white.rds")

