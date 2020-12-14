# 1.0 Load libraries ----

installed.packages("pacman") 
library("pacman")

pacman::p_load(tidyverse,data.table, tictoc, furrr, future,ggpubr, plotly,
               caret, # partition data
               rcompanion, # calculating R2
               pROC,  # calculating and plot AUC 
               plotROC,
               biglasso)  # calculating and plot AUC as ggplot   

# 2.0 Importing Files ----

dir_diabetes_prs_scores <- "/home/opc/storage/uk_biobank/phenotype_data/prs/t2d"
no_cores <- availableCores() - 1

plan(multisession, workers = no_cores)

diabetes_prs_scores <- list.files(dir_diabetes_prs_scores , full.names = T) %>%
   future_map(fread) %>%
   reduce(full_join, by= c("FID" = "FID")) %>%
   select(matches("FID"), matches("Pt_1"))

# # 3.0 Calculating score sum ----
# 
diabetes_prs <-  diabetes_prs_scores %>%
   mutate( prs_diabetes = rowSums(diabetes_prs_scores %>% dplyr::select(- FID), na.rm = TRUE)) %>%
   dplyr::select( FID,prs_diabetes) %>%
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

prev_diabetes_scores <- prs_pheno(diabetes_final,diabetes_prs)

predict_auc <- function(data){
   set.seed(123)
   inTrain<- createDataPartition(y=data$pheno, p=0.75, list=FALSE)
   training<-data[inTrain,]

   testing<-data[-inTrain,]

   glm(pheno ~
          prs_diabetes +
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

predict_auc(prev_diabetes_scores)

tic()
feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
                          "all_chr_genepy_score.feather")
exome_scores <- arrow::read_feather(feather_file) %>% as.data.frame()
toc()

conversion_table <- fread("/home/opc/storage/output/output/Exome_scores/conversion_tables")

exome_scores %>%
   inner_join(conversion_table %>% dplyr::select(V1,V7), by = c('IID' = 'V1')) %>%
   inner_join(prev_diabetes_scores %>% select(eid, pheno
                                         #      ,
                                         # prs_diabetes, sex_f31_0_0,
                                         # age_when_attended_assessment_centre_f21003_0_0,
                                         # genetic_principal_components_f22009_0_1,
                                         # genetic_principal_components_f22009_0_2,
                                         # genetic_principal_components_f22009_0_3,
                                         # genetic_principal_components_f22009_0_4
                                         )
              , by = c('V7' = 'eid')) %>%
   mutate(pheno = ifelse(pheno == 1, "yes", "no")) -> prs_gene_score_diabetes

prs_gene_score_diabetes <- prs_gene_score_diabetes %>%
   select(IID,pheno,V7,everything())

tic()
laso_select <- function(data){
#create training data
    set.seed(1)
    training <- createDataPartition(prs_gene_score_diabetes$pheno, p = 0.75, list=FALSE)

    trainData <- prs_gene_score_diabetes[training,]
#preprocessing
    removeZeroVar3 <- function(df){
       df[, !sapply(df, function(x) min(x) == max(x))]
    }

    print("remove zero variance")
    removeZeroVar3(trainData) -> trainData_fltrd

    print("scale data")
    trainData_fltrd[, -c(1,2,3)] <- lapply( trainData_fltrd[, -c(1,2,3)], function(x) c(scale(x)))

    preprocessed_data <- trainData_fltrd %>%
       mutate(pheno = ifelse(pheno == "yes" , 1, 0),
              pheno = as.factor(as.numeric(pheno)))

    feather_file <- file.path("/home/opc/storage/output/output/Exome_scores/preprocessed_diabetes.feather")
    arrow::write_feather(preprocessed_data, feather_file)

    train.bm <- as.big.matrix(preprocessed_data %>% select(-IID,-pheno,-V7))



   fit.bl1 <- cv.biglasso(train.bm, preprocessed_data$pheno, penalty="lasso", eval.metric="default",
                          family="binomial" ,trace = TRUE, screen = "None",
                          seed = 1234, ncores = 10, nfolds = 10)

}

tic()
feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
                          "preprocessed_diabetes.feather")
preprocessed_diabetes<- arrow::read_feather(feather_file) %>% as.data.frame()
train.bm <- as.big.matrix(preprocessed_diabetes %>% select(-IID,-pheno,-V7))
toc()
tic()
laso_select <- function(data){
   fit.bl1 <- cv.biglasso(train.bm, preprocessed_diabetes$pheno, penalty="lasso", eval.metric="default",
                          family="binomial" ,trace = TRUE, screen = "None",
                          seed = 1234, ncores = 10, nfolds = 10)
}
toc()

tic()

laso_select(preprocessed_diabetes) -> laso_model_diabetes
toc()

geneImp <- function(object, lambda = NULL, ...) {
   
   beta <- predict(object, s = lambda, type = "coef")
   if(is.list(beta)) {
      out <- do.call("cbind", lapply(beta, function(x) x[,1]))
      out <- as.data.frame(out, stringsAsFactors = TRUE)
   } else out <- data.frame(Overall = beta[,1])
   out <- abs(out[rownames(out) != "(Intercept)",,drop = FALSE])
   out
}

geneImp(laso_model_diabetes, lambda = laso_model_diabetes$lambda.min ,scale = TRUE) %>% 
   as.data.frame() %>% 
   rownames_to_column() %>%
   arrange(desc(Overall)) %>% 
   filter(Overall != 0) %>% 
   select(1) -> varbs2
varbs2$rowname -> varbs2


print("saving files")

saveRDS(varbs2, file.path("/home/opc/storage/output/output/prs_genepy/diabetes_genes.rds"))

varbs3 <- coef(laso_model_diabetes)[which(coef(laso_model_diabetes) != 0)]

saveRDS(varbs3, file.path("/home/opc/storage/output/output/prs_genepy/diabetes_genes2.rds"))

varbs2 <- readRDS("/home/opc/storage/output/output/prs_genepy/diabetes_genes.rds")

 score_interest <- prs_gene_score_diabetes %>%
    dplyr::select(varbs2,pheno,IID,V7)
 score_interest <- score_interest %>%
       mutate( exome_score_sum =
               rowSums(scale(score_interest  %>% dplyr::select(-pheno,-IID,-V7)))) %>%
    left_join(prev_diabetes_scores %>%
                 select(eid,sex_f31_0_0,
                        age_when_attended_assessment_centre_f21003_0_0,
                        genetic_principal_components_f22009_0_1,
                        genetic_principal_components_f22009_0_2,
                        genetic_principal_components_f22009_0_3,
                        genetic_principal_components_f22009_0_4,
                        prs_diabetes),  by = c('V7' = 'eid')) %>%
    mutate(pheno = as.factor(ifelse(pheno == "yes" , 1, 0)))


 score_interest <- score_interest %>%
    mutate( exome_score_sum =
               rowSums(scale(score_interest %>%
                          dplyr::select(-pheno,-IID,-sex_f31_0_0,
                                        -age_when_attended_assessment_centre_f21003_0_0,
                                        -genetic_principal_components_f22009_0_1,
                                        -genetic_principal_components_f22009_0_2,
                                        -genetic_principal_components_f22009_0_3,
                                        -genetic_principal_components_f22009_0_4,
                                        -prs_diabetes)),
                                      na.rm = TRUE))
 set.seed(010)
 training <- createDataPartition(score_interest$pheno, p = 0.75, list=FALSE)

 trainData <- score_interest[training,]
 testData <- score_interest[-training,]

 
 glm(as.factor(pheno) ~  prs_diabetes + exome_score_sum +
        age_when_attended_assessment_centre_f21003_0_0
     +sex_f31_0_0
     + genetic_principal_components_f22009_0_1
     + genetic_principal_components_f22009_0_2
     +genetic_principal_components_f22009_0_3
     +genetic_principal_components_f22009_0_4,

     data=trainData, family = binomial) -> model_1



 glm(as.factor(pheno) ~  prs_diabetes
     + age_when_attended_assessment_centre_f21003_0_0
     +sex_f31_0_0
     + genetic_principal_components_f22009_0_1
     + genetic_principal_components_f22009_0_2
     +genetic_principal_components_f22009_0_3
     +genetic_principal_components_f22009_0_4,
     data=trainData, family = binomial) -> model_2

 glm(as.factor(pheno) ~ exome_score_sum +
        age_when_attended_assessment_centre_f21003_0_0
     +sex_f31_0_0
     + genetic_principal_components_f22009_0_1
     + genetic_principal_components_f22009_0_2
     +genetic_principal_components_f22009_0_3
     +genetic_principal_components_f22009_0_4,
     data=trainData, family = "binomial") -> model_3



 predict(model_1, testData) -> prob_1
 predict(model_2, testData) -> prob_2
 predict(model_3, testData) -> prob_3


 testData$prob_1=prob_1
 testData$prob_2=prob_2
 testData$prob_3=prob_3


 g_1 <- roc(pheno ~ prob_1, data = testData, plot=TRUE,
          print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

 g_2 <- roc(pheno ~ prob_2, data = testData, plot=TRUE,
            print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

 g_3 <- roc(pheno ~ prob_3, data = testData, plot=TRUE,
            print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
            
            