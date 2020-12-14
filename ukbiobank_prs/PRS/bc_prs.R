pacman::p_load(tidyverse,data.table, tictoc, furrr, future,ggpubr, plotly,
               caret, # partition data
               rcompanion, # calculating R2
               pROC,  # calculating and plot AUC 
               plotROC,
               doMC,
               glmnet,
               biglasso)  # calculating and plot AUC as ggplot 
dir_scores_khera <- "/home/opc/storage/output/output/score_khera"

no_cores <- 5
plan(multisession, workers = no_cores)

scores_khera <- list.files(dir_scores_khera,pattern = "sscore$", full.names = T) %>%
   future_map(fread) %>%
   reduce(full_join, by= c("#FID" = "#FID")) %>%
   select(matches("#FID"), matches("SCORE"))

scores_sum_khrea <-  scores_khera %>%
   mutate( sum =  rowSums(scores_khera %>% dplyr::select(- `#FID`), na.rm = TRUE)) %>%
   dplyr::select(`#FID`,sum) %>%
   rename( id = `#FID`)


prev_scores_sum_khera <- prs_pheno(scores_sum_khrea, "prevalent")

# tic()
# feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
#                           "preprocessed_traindata.feather")
# preprocessed_data <- arrow::read_feather(feather_file) %>% as.data.frame()
# toc()
#
# bc_final <- preprocessed_data  %>%
#                 inner_join(prev_scores_sum_khera %>% select(id, pheno),
#                             by = c('V7' = 'id')) %>%
#             mutate(pheno = as.factor(pheno))

# tic()
# feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
#                           "preprocessed_bc.feather")
# arrow::write_feather(bc_final,feather_file)
# toc()

tic()
feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
                          "preprocessed_bc.feather")
preprocessed_bc<- arrow::read_feather(feather_file) %>% as.data.frame()
toc()


train.bm <- as.big.matrix(preprocessed_bc %>% select(-IID,-pheno,-V7))
toc()

tic()
laso_select <- function(data){
   laso_model_bc <- cv.biglasso(data, preprocessed_bc$pheno, penalty="lasso", 
                                eval.metric="default",family="binomial" ,
                                trace = TRUE, screen = "None",
                                seed = 12, ncores = 10, nfolds = 10)
}
toc()

tic()

laso_select(train.bm) -> laso_model_bc

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

geneImp(laso_model_bc, lambda = laso_model_bc$lambda.min ,scale = TRUE) %>% 
   as.data.frame() %>% 
   rownames_to_column() %>%
   arrange(desc(Overall)) %>% 
   filter(Overall != 0) %>% 
   select(1) -> varbs2
varbs2$rowname -> varbs2


print("saving files")

saveRDS(varbs2, file.path("/home/opc/storage/output/output/prs_genepy/bc_genes2.rds"))

readRDS(file.path("/home/opc/storage/output/output/prs_genepy/bc_genes2.rds")) -> bc_genes

tic()
feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
                          "all_chr_genepy_score.feather")
exome_scores <- arrow::read_feather(feather_file) %>% as.data.frame()
toc()
# 
conversion_table <- fread("/home/opc/storage/output/output/Exome_scores/conversion_tables")

exome_scores %>%
   inner_join(conversion_table %>% dplyr::select(V1,V7), by = c('IID' = 'V1')) %>%
   inner_join(prev_scores_sum_khera %>% select(id, pheno
                                         #      ,
                                         # prs_diabetes, sex_f31_0_0,
                                         # age_when_attended_assessment_centre_f21003_0_0,
                                         # genetic_principal_components_f22009_0_1,
                                         # genetic_principal_components_f22009_0_2,
                                         # genetic_principal_components_f22009_0_3,
                                         # genetic_principal_components_f22009_0_4
   )
   , by = c('V7' = 'id')) %>%
   mutate(pheno = ifelse(pheno == 1, "yes", "no")) -> prs_gene_score_bc

score_interest <- prs_gene_score_bc %>%
   dplyr::select(bc_genes,pheno,IID,V7)

score_interest <- score_interest %>%
   mutate( exome_score_sum =
              rowSums(scale(score_interest  %>% dplyr::select(-pheno,-IID,-V7)))) %>%
   left_join(prev_scores_sum_khera %>%
                select(id,
                       age,
                       genetic_principal_components_f22009_0_1,
                       genetic_principal_components_f22009_0_2,
                       sum),  by = c('V7' = 'id')) %>%
   mutate(pheno = as.factor(ifelse(pheno == "yes" , 1, 0)))
#
#
set.seed(12)
training <- createDataPartition(score_interest$pheno, p = 0.75, list=FALSE)

cad_trainData <- score_interest[training,]
cad_testData <- score_interest[-training,]

glm(as.factor(pheno) ~  sum + exome_score_sum +
       age

    + genetic_principal_components_f22009_0_1
    + genetic_principal_components_f22009_0_2,

    data=cad_trainData, family = binomial) -> model_1

glm(as.factor(pheno) ~  sum
    + age
    + genetic_principal_components_f22009_0_1
    + genetic_principal_components_f22009_0_2,
    data=cad_trainData, family = binomial) -> model_2

glm(as.factor(pheno) ~ exome_score_sum +
       age
    + genetic_principal_components_f22009_0_1
    + genetic_principal_components_f22009_0_2,
    data=cad_trainData, family = binomial) -> model_3



predict(model_1, cad_testData) -> prob_1
predict(model_2, cad_testData) -> prob_2
predict(model_3, cad_testData) -> prob_3


cad_testData$prob_1 = prob_1
cad_testData$prob_2 = prob_2
cad_testData$prob_3 = prob_3


g_1 <- roc(pheno ~ prob_1, data = cad_testData, plot=TRUE,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

g_2 <- roc(pheno ~ prob_2, data = cad_testData, plot=TRUE,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

g_3 <- roc(pheno ~ prob_3, data = cad_testData, plot=TRUE,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

longtest <- melt_roc(cad_testData, "pheno", c("prob_1", "prob_2","prob_3")) %>%
   mutate(name = case_when(
      name == "prob_1" ~ "gps + prs",
      name == "prob_2" ~"prs",
      name == "prob_3" ~"gps"))




model <- unique(longtest$name)
model_info <- data.frame(model,
                         group = rank(model))

ggplot(longtest, aes( m = M, d = as.factor(as.numeric(D)), color = name)) +
   geom_roc(n.cuts=20,labels=FALSE) -> p

left_join(model_info, calc_auc(p)) %>%
   select(-group, -PANEL) %>%
   arrange(desc(AUC)) %>%
   mutate( AUC = round(AUC,3)) -> aucs



p +
   annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()),
                     xmin=0.5, xmax=1, ymin=0, ymax=0.5) +
   theme_pubr()


         glm(pheno ~ sum
             + genetic_principal_components_f22009_0_1
             + genetic_principal_components_f22009_0_2,
             data=cad_trainData, family = "binomial") -> model_4

         glm(pheno ~ genetic_principal_components_f22009_0_1
             + genetic_principal_components_f22009_0_2,
             data=cad_trainData, family = "binomial") -> model_5

         predict(model_4, cad_testData) -> prob_4
         predict(model_5, cad_testData) -> prob_5

         cad_testData$prob_4=prob_4
         cad_testData$prob_5=prob_5

         g_4 <- roc(pheno ~ prob_4, data = cad_testData, plot=TRUE,
                    print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

         g_5 <- roc(pheno ~ prob_5, data = cad_testData, plot=TRUE,
                    print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)
         