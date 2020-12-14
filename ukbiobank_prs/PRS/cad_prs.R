# # 1.0 Load libraries ----
# 
# installed.packages("pacman") 
# library("pacman")
# 
pacman::p_load(tidyverse,data.table, tictoc, furrr, future,ggpubr, plotly,
               caret, # partition data
               rcompanion, # calculating R2
               pROC,  # calculating and plot AUC 
               plotROC,
               doMC,
               glmnet,
               biglasso)  # calculating and plot AUC as ggplot   

# # 2.0 Importing Files ----

# dir_cad_prs_scores <- "/home/opc/storage/uk_biobank/phenotype_data/prs/cad"
# no_cores <- 5
# 
# plan(multisession, workers = no_cores)
# 
# cad_prs_scores <- list.files(dir_cad_prs_scores , full.names = T) %>%
#    future_map(fread) %>%
#    reduce(full_join, by= c("FID" = "FID")) %>%
#    select(matches("FID"), matches("Pt_1"))
# 
# # # 3.0 Calculating score sum ----
# #
# cad_prs <-  cad_prs_scores %>%
#    mutate( prs_cad = rowSums(cad_prs_scores %>% dplyr::select(- FID), na.rm = TRUE)) %>%
#    dplyr::select( FID,prs_cad) %>%
#    rename( eid = `FID`)
# 
# # #
# # # ## function
# # #
# prs_pheno <- function(data1, data2){
# 
#    cases_ids <- data1 %>%
#       filter( prev_inc == "prevalent") %>%
#       dplyr::select(eid,age,sex_f31_0_0,age_when_attended_assessment_centre_f21003_0_0,
#                     genetic_principal_components_f22009_0_1,
#                     genetic_principal_components_f22009_0_2,
#                     genetic_principal_components_f22009_0_3,
#                     genetic_principal_components_f22009_0_4)  %>%
#       distinct(.keep_all = TRUE) %>%
#       mutate(pheno = 1)
# 
#    control_ids <-
#       ukbb_df %>% dplyr::select(eid,sex_f31_0_0,ethnic_background_f21000_0_0,
#                                 outliers_for_heterozygosity_or_missing_rate_f22027_0_0,
#                                 sex_chromosome_aneuploidy_f22019_0_0,
#                                 genetic_sex_f22001_0_0,
#                                 age_when_attended_assessment_centre_f21003_0_0,
#                                 genetic_principal_components_f22009_0_1,
#                                 genetic_principal_components_f22009_0_2,
#                                 genetic_principal_components_f22009_0_3,
#                                 genetic_principal_components_f22009_0_4) %>%
#       filter(ethnic_background_f21000_0_0 == "British") %>%
#       filter(outliers_for_heterozygosity_or_missing_rate_f22027_0_0 != "Yes" |
#                 is.na(outliers_for_heterozygosity_or_missing_rate_f22027_0_0)) %>%
#       filter(sex_chromosome_aneuploidy_f22019_0_0 != "Yes" |
#                 is.na(sex_chromosome_aneuploidy_f22019_0_0)) %>%
#       mutate( sex_check= case_when(
#          sex_f31_0_0 == genetic_sex_f22001_0_0 ~ TRUE,
#          TRUE ~ FALSE
#       )) %>%
#       filter(sex_check == TRUE) %>%
#       dplyr::select(-sex_check) %>%
#       anti_join(cases_ids) %>%
#       dplyr::select(eid,age_when_attended_assessment_centre_f21003_0_0,sex_f31_0_0,
#                     genetic_principal_components_f22009_0_1,
#                     genetic_principal_components_f22009_0_2,
#                     genetic_principal_components_f22009_0_3,
#                     genetic_principal_components_f22009_0_4)  %>%
#       mutate(age = age_when_attended_assessment_centre_f21003_0_0) %>%
#       distinct(.keep_all = TRUE) %>%
#       mutate(pheno = 0)
# 
#    all_ids <- cases_ids %>%
#       bind_rows(control_ids)
# 
# 
# #    #ukbb_df %>%
# #    #   dplyr::select(eid, matches("illnesses_of_mother")) %>%
# #    #   gather("relationship" , "illness", - eid) %>%
# #    #   filter(str_detect(illness, 'Breast')) %>%
# #    #   dplyr::select(-relationship) %>%
# #    #   distinct(eid, .keep_all = TRUE)-> fam_hist
# #
# scores_sum_pheno <- data2 %>%
#    left_join(all_ids, by = c("eid"="eid")) %>%
#    filter(!is.na(pheno)) #%>%
#          # left_join(fam_hist, by = c("id"="eid")) %>%
#          # mutate(ill = ifelse(illness == "Breast cancer", 1, 0 )) %>%
#          # mutate(ill = replace_na(ill,0)) %>%
#          # left_join(ukbb_pc_df %>% dplyr::select(eid, genetic_principal_components_f22009_0_1,
#          #                                       genetic_principal_components_f22009_0_2  ),
#          #          by = c("id" ="eid"))
#    return(scores_sum_pheno)
# }
# 
# prev_cad_scores <- prs_pheno(cad_final,cad_prs)
# 
# predict_auc <- function(data){
#    set.seed(123)
#    inTrain<- createDataPartition(y=data$pheno, p=0.75, list=FALSE)
#    training<-data[inTrain,]
# 
#    testing<-data[-inTrain,]
# 
#    glm(pheno ~
#           prs_cad +
#           age_when_attended_assessment_centre_f21003_0_0 +
#           sex_f31_0_0 +
#           genetic_principal_components_f22009_0_1 +
#           genetic_principal_components_f22009_0_2 +
#           genetic_principal_components_f22009_0_3 +
#           genetic_principal_components_f22009_0_4 ,
#        training , family = binomial("logit") , maxit = 100) -> model
# 
#    predict(model, testing) -> prob
#    testing$prob=prob
#    # calculating AUC
#    g <- roc(pheno ~ prob, data = testing, plot=TRUE,
#             print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
# }
# 
# predict_auc(prev_cad_scores)
# #
# # 
# 
 tic()
feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
                          "preprocessed_internal_maf_cad.feather")
preprocessed_cad <- arrow::read_feather(feather_file) %>% as.data.frame()
toc()
# # # # 
# 
# tic()
# feather_file <- file.path("/home/opc/storage/output/output/Exome_scores/200K_internal_MAF",
#                           "200k_exome_wide_scores_internal_MAF.feather")
# exome_scores <- arrow::read_feather(feather_file) %>% as.data.frame()
# toc()

# tic()
# feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
#                           "all_chr_genepy_score.feather")
# exome_scores <- arrow::read_feather(feather_file) %>% as.data.frame()
# toc()
# # 
# conversion_table <- fread("/home/opc/storage/output/output/Exome_scores/conversion_tables")
# 
# exome_scores %>%
#    inner_join(conversion_table %>% dplyr::select(V1,V7), by = c('IID' = 'V1')) %>%
#    inner_join(prev_cad_scores %>% select(eid, pheno
#                                          #      ,
#                                          # prs_diabetes, sex_f31_0_0,
#                                          # age_when_attended_assessment_centre_f21003_0_0,
#                                          # genetic_principal_components_f22009_0_1,
#                                          # genetic_principal_components_f22009_0_2,
#                                          # genetic_principal_components_f22009_0_3,
#                                          # genetic_principal_components_f22009_0_4
#                                          )
#               , by = c('V7' = 'eid')) %>%
#    mutate(pheno = ifelse(pheno == 1, "yes", "no")) -> prs_gene_score_cad
# 
# prs_gene_score_cad <- prs_gene_score_cad %>%
#    select(IID,pheno,V7,everything())
# 
# tic()
# laso_select <- function(data){
# #create training data
#     set.seed(1)
#     training <- createDataPartition(prs_gene_score_cad$pheno, p = 0.75, list=FALSE)
# 
#     trainData <- prs_gene_score_cad[training,]
# #preprocessing
#     removeZeroVar3 <- function(df){
#        df[, !sapply(df, function(x) min(x) == max(x))]
#     }
# 
#     print("remove zero variance")
#     removeZeroVar3(trainData) -> trainData_fltrd
#     trainData_fltrd <- trainData_fltrd %>%
#         select(IID,pheno,V7,everything())
#     print("scale data")
#     trainData_fltrd[, -c(1,2,3)] <- lapply( trainData_fltrd[, -c(1,2,3)], function(x) c(scale(x)))
# 
#     preprocessed_data <- trainData_fltrd %>%
#        mutate(pheno = ifelse(pheno == "yes" , 1, 0),
#               pheno = as.factor(as.numeric(pheno)))
# 
# #     preprocessed_data <- trainData_fltrd %>% 
# #        mutate(pheno = ifelse(pheno == "yes" , 1, 0),
# #               pheno = as.factor(as.numeric(pheno)))
# # 
     # feather_file <- file.path("/home/opc/storage/output/output/Exome_scores/preprocessed_internal_maf_cad.feather")
     # arrow::write_feather(preprocessed_data, feather_file)
#toc()
#
#
# genepy_exome_scores %>%
#  inner_join(conversion_table %>% dplyr::select(V1,V7), by = c('IID' = 'V1')) %>%
#  inner_join(prev_cad_scores %>% select(eid, pheno, prs_cad, sex_f31_0_0,
#                                        age_when_attended_assessment_centre_f21003_0_0,
#                                        genetic_principal_components_f22009_0_1,
#                                        genetic_principal_components_f22009_0_2,
#                                        genetic_principal_components_f22009_0_3,
#                                        genetic_principal_components_f22009_0_4)
#                , by = c('V7' = 'eid')) %>%
#     mutate(pheno = ifelse(pheno == 1, "yes", "no")) -> prs_gene_score_cad

# lasso function for gene selection
# print("start lasso")
# 
# laso_select <- function(data){
#    
#    
#    #testData <- data[-training,]
#    #feather_file <- file.path("/home/opc/storage/output/output/prs_genepy/preprocessed_cad.feather")
#    #arrow::write_feather(pre_proc_trainData, feather_file)
#    
#    fitControl <- trainControl(method = "repeatedcv",
#                               number = 2,
#                               repeats = 2,
#                               classProbs = TRUE,
#                               summaryFunction = twoClassSummary,
#                               returnData = FALSE,
#                               returnResamp = "final",
#                               trim=TRUE
#    )
#    
#    set.seed(1)
#    myGrid <- expand.grid(
#       alpha = 1,
#       lambda = 10^seq(1, -7, length = 10))
#    
#    # library(doParallel)
#    # cl <- makePSOCKcluster(2)
#    # registerDoParallel(cl)
#    print("start computing lasso")
#    lassoModel20 <- train(x=data %>% dplyr::select( -IID, -pheno),
#                          y=data$pheno,
#                          method = "glmnet",  metric="ROC",
#                          trControl = fitControl,
#                          verbose=TRUE, tuneGrid = myGrid)
# }
#    #,
#    #                         preProcess = c("zv","medianImpute","center","scale"),
#    #                       model = FALSE)
# 
# 
# 
# tic()
# laso_select(preprocessed_cad) -> laso_model_cad
# toc()
# varImp(laso_model_cad, scale = TRUE)$importance %>% 
#    as.data.frame() %>% 
#    rownames_to_column() %>%
#    arrange(desc(Overall)) %>% 
#    filter(Overall != 0) %>% 
#    select(1) -> varbs2
# varbs2$rowname -> varbs2
# 
# saveRDS(varbs2, file.path("/home/opc/storage/output/output/prs_genepy/cad_genes.rds"))
# ## plot gps+prs
# 
# 
#    glm(pheno ~  prs_cad + exome_score_sum
#        + age_when_attended_assessment_centre_f21003_0_0
#        +sex_f31_0_0
#        + genetic_principal_components_f22009_0_1
#        + genetic_principal_components_f22009_0_2
#        +genetic_principal_components_f22009_0_3
#        +genetic_principal_components_f22009_0_4,
# 
#        data=training_cad, family = binomial) -> model_1
# 
#    #model_1 %>% stepAIC(direction='both',trace = FALSE) -> model_2
# 
#    glm(pheno ~ prs_cad + age_when_attended_assessment_centre_f21003_0_0
#        +sex_f31_0_0
#        + genetic_principal_components_f22009_0_1
#        + genetic_principal_components_f22009_0_2
#        +genetic_principal_components_f22009_0_3
#        +genetic_principal_components_f22009_0_4,
#        data=training_cad, family = binomial) -> model_2
# 
#    glm(pheno ~  exome_score_sum + age_when_attended_assessment_centre_f21003_0_0
#        +sex_f31_0_0
#        + genetic_principal_components_f22009_0_1
#        + genetic_principal_components_f22009_0_2
#        +genetic_principal_components_f22009_0_3
#        +genetic_principal_components_f22009_0_4,
#        data=training_cad, family = "binomial") -> model_3
# 
# 
# 
#    predict(model_1, testing_cad) -> prob_1
#    predict(model_2, testing_cad) -> prob_2
#    predict(model_3, testing_cad) -> prob_3
# 
# 
#    testing_cad$prob_1=prob_1
#    testing_cad$prob_2=prob_2
#    testing_cad$prob_3=prob_3
# 
#    g <- roc(pheno ~ prob_1, data = testing_cad, plot=TRUE,
#             print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
# 
#    longtest <- melt_roc(testing_cad, "pheno", c("prob_1", "prob_2","prob_3")) %>%
#       mutate(name = case_when(
#          name == "prob_1" ~ "gps + prs",
#          name == "prob_2" ~"prs",
#          name == "prob_3" ~"gps"))
# 
#    model <- unique(longtest$name)
#    model_info <- data.frame(model,
#                             group = rank(model))
# 
#    ggplot(longtest, aes(d = D, m = M, color = name)) +
#       geom_roc(n.cuts=20,labels=FALSE) ->p
# 
#    left_join(model_info, calc_auc(p)) %>%
#       select(-group, -PANEL) %>%
#       arrange(desc(AUC)) %>%
#       mutate( AUC = round(AUC,3)) -> aucs
# 
# 
# 
#    p +
#       annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()),
#                         xmin=0.5, xmax=1, ymin=0, ymax=0.5) +
#       theme_pubr()
# 
# tic()  
# pre_proc_values <- preProcess(trainData, method = c("zv","medianImpute","center","scale"))
# toc()
# tic()
# pre_proc_trainData <- predict(pre_proc_values, trainData)
# toc()



# nearZeroVar(samp_train, saveMetrics = TRUE, foreach = TRUE,
#             allowParallel = TRUE)  -> emad
# 
# samp_train_fltrd <- samp_train[, emad$zeroVar == 'FALSE']
# 
# removeZeroVar3 <- function(df){
#    df[, !sapply(df, function(x) min(x) == max(x))]
# }
# # 
# samp_train <- trainData %>% tail(100)
# removeZeroVar3(samp_train) -> samp_train_fltrd
# 
# samp_train_fltrd[, -c(1,2)] <- lapply( samp_train_fltrd[, -c(1,2)], function(x) c(scale(x)))
# #samp_train_fltrd[, -c(1,2)] <- scale(samp_train_fltrd[, -c(1,2)])
#samp_train_fltrd[, -c(1)] <- scale(samp_train_fltrd[, -c(1)])
# tic()
# pre_proc_values <- preProcess(samp_train, method = c("zv","scale"))
# toc()
# print("Apply preprocess to the data sets")
# tic()
# pre_proc_trainData <- predict(pre_proc_values, samp_train)
# toc()

train.bm <- as.big.matrix(preprocessed_cad %>% select(-IID,-pheno,-V7))


tic()
laso_select <- function(data){
fit.bl1 <- cv.biglasso(data, preprocessed_cad$pheno, penalty="lasso", eval.metric="default",
                      family="binomial" ,trace = TRUE, screen = "None",
                      seed = 1234, ncores = 10, nfolds = 10)
}
toc()

# laso_select <- function(data){
#    registerDoMC(cores = 10)
# fit.bl <- cv.glmnet(data.matrix(data %>% select(-IID,-pheno)), data$pheno,  alpha =1,
#                     type.measure = "auc",  family="binomial",trace.it = 1,
#                     parallel = TRUE)
# }

print("starts lasso")

tic()
laso_select(train.bm) -> cad_model
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


geneImp(cad_model, lambda = cad_model$lambda.min ,scale = TRUE) %>%
   as.data.frame() %>%
   rownames_to_column() %>%
   arrange(desc(Overall)) %>%
   filter(Overall != 0) %>%
   select(1) -> varbs2
varbs2$rowname -> varbs2

print("saving files")

saveRDS(varbs2, file.path("/home/opc/storage/output/output/prs_genepy/cad_genes_internal_maf_cad.rds"))

# # formula <- as.formula(pheno ~ .)
# # X <- model.matrix(formula, samp_train)
# 
# # tic()
# # feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
# #                           "prs_gene_score_cad.feather")
# # prs_cad<- arrow::read_feather(feather_file) %>% as.data.frame()Z
# # toc()
# # 
# readRDS(file.path("/home/opc/storage/output/output/prs_genepy/cad50_genes2.rds")) -> cad_genes
# 
# score_interest <- prs_cad %>%
#    dplyr::select(cad_genes,pheno,IID,V7)
# # 
# score_interest <- score_interest %>%
#    mutate( exome_score_sum =
#               rowSums(scale(score_interest  %>% dplyr::select(-pheno,-IID,-V7)))) %>%
#    left_join(prev_cad_scores %>%
#                 select(eid,sex_f31_0_0,
#                        age_when_attended_assessment_centre_f21003_0_0,
#                        genetic_principal_components_f22009_0_1,
#                        genetic_principal_components_f22009_0_2,
#                        genetic_principal_components_f22009_0_3,
#                        genetic_principal_components_f22009_0_4,
#                        prs_cad),  by = c('V7' = 'eid')) %>%
#    mutate(pheno = as.factor(ifelse(pheno == "yes" , 1, 0)))
# # # 
# # # 
# set.seed(1)
# training <- createDataPartition(score_interest$pheno, p = 0.75, list=FALSE)
# 
# cad_trainData <- score_interest[training,]
# cad_testData <- score_interest[-training,]
# 
# glm(as.factor(pheno) ~  prs_cad + exome_score_sum +
#        age_when_attended_assessment_centre_f21003_0_0
#     +sex_f31_0_0
#     + genetic_principal_components_f22009_0_1
#     + genetic_principal_components_f22009_0_2
#     +genetic_principal_components_f22009_0_3
#     +genetic_principal_components_f22009_0_4,
# 
#     data=cad_trainData, family = binomial) -> model_1
# # 
# glm(as.factor(pheno) ~  prs_cad
#     + age_when_attended_assessment_centre_f21003_0_0
#     + sex_f31_0_0
#     + genetic_principal_components_f22009_0_1
#     + genetic_principal_components_f22009_0_2
#     + genetic_principal_components_f22009_0_3
#     + genetic_principal_components_f22009_0_4,
#     data=cad_trainData, family = binomial) -> model_2
# 
# glm(as.factor(pheno) ~ exome_score_sum +
#        age_when_attended_assessment_centre_f21003_0_0
#     + sex_f31_0_0
#     + genetic_principal_components_f22009_0_1
#     + genetic_principal_components_f22009_0_2
#     + genetic_principal_components_f22009_0_3
#     + genetic_principal_components_f22009_0_4,
#     data=cad_trainData, family = binomial) -> model_3
# 
# 
# 
# predict(model_1, cad_testData) -> prob_1
# predict(model_2, cad_testData) -> prob_2
# predict(model_3, cad_testData) -> prob_3
# 
# 
# cad_testData$prob_1 = prob_1
# cad_testData$prob_2 = prob_2
# cad_testData$prob_3 = prob_3
# 
# 
# g_1 <- roc(pheno ~ prob_1, data = cad_testData, plot=TRUE,
#            print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
# 
# g_2 <- roc(pheno ~ prob_2, data = cad_testData, plot=TRUE,
#            print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
# 
# g_3 <- roc(pheno ~ prob_3, data = cad_testData, plot=TRUE,
#            print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
# 
#    longtest <- melt_roc(cad_testData, "pheno", c("prob_1", "prob_2","prob_3")) %>%
#       mutate(name = case_when(
#          name == "prob_1" ~ "gps + prs",
#          name == "prob_2" ~"prs",
#          name == "prob_3" ~"gps"))
# 
# 
# 
# 
#    model <- unique(longtest$name)
#    model_info <- data.frame(model,
#                             group = rank(model))
# 
#    ggplot(longtest, aes( m = M, d = as.factor(as.numeric(D)), color = name)) +
#       geom_roc(n.cuts=20,labels=FALSE) -> p
# 
#    left_join(model_info, calc_auc(p)) %>%
#       select(-group, -PANEL) %>%
#       arrange(desc(AUC)) %>%
#       mutate( AUC = round(AUC,3)) -> aucs
# 
# 
# 
#    p +
#       annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()),
#                         xmin=0.5, xmax=1, ymin=0, ymax=0.5) +
#       theme_pubr()
# 
# cad_gps_prs <- function(training,testing){
#       glm(pheno ~ prs_cad
#           + genetic_principal_components_f22009_0_1
#           + genetic_principal_components_f22009_0_2
#           + genetic_principal_components_f22009_0_3
#           + genetic_principal_components_f22009_0_4,
#           data=training, family = "binomial") -> model_4
# 
#       glm(pheno ~ exome_score_sum
#           + genetic_principal_components_f22009_0_1
#           + genetic_principal_components_f22009_0_2
#           + genetic_principal_components_f22009_0_3
#           + genetic_principal_components_f22009_0_4,
#           data=training, family = "binomial") -> model_5
# 
#       predict(model_4, testing) -> prob_4
#       predict(model_5, testing) -> prob_5
# 
#       testing$prob_4=prob_4
#       testing$prob_5=prob_5
# 
#       g_4 <- roc(pheno ~ prob_4, data = cad_testData, plot=TRUE,
#                  print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
# 
#       g_5 <- roc(pheno ~ prob_5, data = cad_testData, plot=TRUE,
#                  print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
# 
#       longtest2 <- melt_roc(testing, "pheno", c("prob_4", "prob_5")) %>%
#          mutate(name = case_when(
#             name == "prob_4" ~ "prs only",
#             name == "prob_5" ~"gps only"))
# 
#       model <- unique(longtest2$name)
#       model_info <- data.frame(model,
#                                group = rank(model))
# 
#       ggplot(longtest2, aes(d = as.factor(as.numeric(D)), m = M, color = name)) +
#          geom_roc(n.cuts=20,labels=FALSE) ->p2
# 
#       left_join(model_info, calc_auc(p2)) %>%
#          select(-group, -PANEL) %>%
#          arrange(desc(AUC)) %>%
#          mutate( AUC = round(AUC,3)) -> aucs
#       p2 +
#          annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()),
#                            xmin=0.5, xmax=1, ymin=0, ymax=0.5) +
#          theme_pubr()
#       roc.test(g_4,g_5)
#    }
# 
# cad_gps_prs(cad_trainData,cad_testData)
