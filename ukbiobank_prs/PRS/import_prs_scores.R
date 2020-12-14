# 1.0 Load libraries ----

installed.packages("pacman") 
library("pacman")

pacman::p_load(tidyverse,data.table, tictoc, furrr, future,ggpubr, plotly,
               caret, # partition data
               rcompanion, # calculating R2
               pROC,  # calculating and plot AUC 
               plotROC)  # calculating and plot AUC as ggplot   

# 2.0 Importing Files ----

dir_scores_313 <- "/home/opc/storage/output/output/score_313"
dir_scores_khera <- "/home/opc/storage/output/output/score_khera"
dir_scores_bcac <- "/home/opc/storage/output/output/score_bcac"


no_cores <- availableCores() - 1

plan(multisession, workers = no_cores)

scores_313 <- list.files(dir_scores_313,pattern = "sscore$", full.names = T) %>%
   future_map(fread) %>%
   reduce(full_join, by= c("#FID" = "#FID")) %>%
   select(matches("#FID"), matches("SCORE"))

plan(multisession, workers = no_cores)

scores_khera <- list.files(dir_scores_khera,pattern = "sscore$", full.names = T) %>%
   future_map(fread) %>%
   reduce(full_join, by= c("#FID" = "#FID")) %>%
   select(matches("#FID"), matches("SCORE"))

scores_bcac <- fread(file.path(dir_scores_bcac,"breast_fem_british_bcac_0.2.best")) 

# 3.0 Calculating score sum ----

scores_sum_313 <-  scores_313 %>%
   mutate( sum =  rowSums(scores_313 %>% dplyr::select(- `#FID`), na.rm = TRUE)) %>%
   dplyr::select(`#FID`,sum) %>%
   rename( id = `#FID`)

scores_sum_khrea <-  scores_khera %>%
   mutate( sum =  rowSums(scores_khera %>% dplyr::select(- `#FID`), na.rm = TRUE)) %>%
   dplyr::select(`#FID`,sum) %>%
   rename( id = `#FID`)


prev_scores_sum_313 <- prs_pheno(scores_sum_313, "prevalent")

inc_scores_sum_313 <- prs_pheno(scores_sum_313, "incident")

all_scores_sum_313 <- prs_pheno_all(scores_sum_313)


prev_scores_sum_khera <- prs_pheno(scores_sum_khrea, "prevalent")

inc_scores_sum_khera <- prs_pheno(scores_sum_khrea, "incident")

all_scores_sum_khera <- prs_pheno_all(scores_sum_khrea)

prev_scores_sum_bcac <- prev_scores_sum_313 %>% select(-sum) %>%  
   inner_join(scores_bcac %>% dplyr::select(FID,PRS), by = c("id" = "FID")) %>% 
   mutate(sum = PRS)
inc_scores_sum_bcac <- inc_scores_sum_313 %>% select(-sum) %>% 
   inner_join(scores_bcac %>% dplyr::select(FID,PRS), by = c("id" = "FID")) %>% 
   mutate(sum = PRS)

all_scores_sum_bcac <-  all_scores_sum_313 %>% select(-sum) %>% 
   inner_join(scores_bcac %>% dplyr::select(FID,PRS), by = c("id" = "FID")) %>% 
   mutate(sum = PRS)

# 4.0 genepy score integration ----

# import files
#fread(file = file.path("/home/opc/storage/output/output/Exome_scores" ,
#                       'model.txt')) -> model5
training<- readRDS(file =  "/home/opc/storage/output/output/prs_genepy/breast_white_training.rds")
testing<- readRDS(file =  "/home/opc/storage/output/output/prs_genepy/breast_white_testing.rds")
training_cad<- readRDS(file =  "/home/opc/storage/output/output/prs_genepy/cad_white_training.rds")
testing_cad<- readRDS(file =  "/home/opc/storage/output/output/prs_genepy/cad_white_testing.rds")
training_diabetes<- readRDS(file =  "/home/opc/storage/output/output/prs_genepy/diabetes_white_training.rds")
testing_diabetes<- readRDS(file =  "/home/opc/storage/output/output/prs_genepy/diabetes_white_testing.rds")
#exome_scores <- fread("/home/opc/storage/output/output/Exome_scores/Genome_Wide_Exome_Score")

#imp_varbs <- readRDS(file =  "/home/opc/storage/output/output/prs_genepy/gps_genes_breast_white.rds")
# tidy files
# model5 <-map(seq(from=1, to=17608, by= 1),
#              possibly(function(x) {
#                 glm(pheno ~ . - IID - V7   , 
#                     
#                     prs_gene_score %>% dplyr::select(x, IID, pheno, V7, genetic_principal_components_f22009_0_1,
#                                                      genetic_principal_components_f22009_0_1, age), 
#                     family = "binomial") %>% 
#                    tidy()}, NA_real_) )  %>% 
#    reduce(rbind) %>%
#    distinct(term, .keep_all = TRUE)
# 
# modell <- model5 %>% 
#    filter(str_starts(term, "\\w+" ), !is.na(p.value), p.value < 0.05)



#exome_scores %>% colnames() %>% make.unique() -> new_col_names
#colnames(exome_scores) <- new_col_names

#conversion_table <- fread("/home/opc/storage/output/output/Exome_scores/conversion_tables")


#exome_scores %>%  
#   inner_join(conversion_table %>% dplyr::select(V1,V7), by = c('IID' = 'V1')) %>% 
#   inner_join(prev_scores_sum_313, by = c('V7' = 'id')) -> prs_gene_score

#exome_scores %>%  
#   inner_join(conversion_table %>% dplyr::select(V1,V7), by = c('IID' = 'V1')) %>% 
#   inner_join(prev_scores_sum_313, by = c('V7' = 'id')) -> prs_gene_score




# prs_gene_score %>% 
#    select(modell$term,pheno,IID,genetic_principal_components_f22009_0_1, 
#           genetic_principal_components_f22009_0_2)  -> score_interest
# 
# # calcultaing AUC 
# 
# set.seed(1)
# inTrain<- createDataPartition(y=score_interest$pheno, p=0.75, list=FALSE)
# training<-score_interest[inTrain,]
# 
# testing<-score_interest[-inTrain,]
# 
# glm(pheno ~  . - IID , data=training, family = "binomial") -> model_1
# glm(pheno ~ sum + age , data=training, family = "binomial") -> model_2
# glm(pheno ~  . - IID - sum ,data=training, family = "binomial") -> model_3
# 
# predict(model_1, testing) -> prob_1
# predict(model_2, testing) -> prob_2
# predict(model_3, testing) -> prob_3
# 
# testing$prob_1=prob_1
# testing$prob_2=prob_2
# testing$prob_3=prob_3
# 
# g_1 <- roc(pheno ~ prob_1, data = testing, plot=TRUE, print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
# g_2 <- roc(pheno ~ prob_2, data = testing, plot=TRUE, print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
# g_3 <- roc(pheno ~ prob_3, data = testing, plot=TRUE,print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
# 
# longtest <- melt_roc(testing, "pheno", c("prob_1", "prob_2","prob_3")) %>% 
#    mutate(name = case_when(
#       name == "prob_1" ~ "genepy + prs", 
#       name == "prob_2" ~"prs",
#       name == "prob_3" ~"genepy"))
# 
# ggplot(longtest, aes(d = D, m = M, color = name)) + 
#    geom_roc(n.cuts=20,labels=FALSE) -> p
# 
# model <- unique(longtest$name)
# model_info <- data.frame(model,
#                          group = rank(model))
# left_join(model_info, calc_auc(p)) %>%
#    select(-group, -PANEL) %>%
#    arrange(desc(AUC)) %>% 
#    mutate( AUC = round(AUC,3)) -> aucs
# 
