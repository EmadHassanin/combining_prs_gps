# 1.0 Load libraries ----

installed.packages("pacman") 
library("pacman")

pacman::p_load(tidyverse,data.table, tictoc, furrr, future,ggpubr, plotly,
               caret, # partition data
               rcompanion, # calculating R2
               pROC,  # calculating and plot AUC 
               plotROC)  # calculating and plot AUC as ggplot   

# 2.0 Importing Files ----

#dir_scores <- "/home/opc/storage/output/output/score_313"
#dir_scores <- "/home/opc/storage/output/output/score_khera"

dir_scores <- "/home/opc/storage/output/output/score_khera"
no_cores <- availableCores() - 1

plan(multisession, workers = no_cores)

scores <- list.files(dir_scores,pattern = "sscore$", full.names = T) %>%
   future_map(fread) %>%
   reduce(full_join, by= c("FID" = "FID")) %>%
   select(matches("#FID"), matches("SCORE"))

# 3.0 Calculating score sum ----

scores_sum <-  scores %>%
   mutate( sum =  rowSums(scale(scores %>% dplyr::select(- `#FID`)), na.rm = TRUE)) %>%
   dplyr::select(`#FID`,sum) %>%
   rename( id = `#FID`)



cases_ids <- breast_fem_white  %>%
   filter( prev_inc == "incident") %>%
   #filter( prev_inc == "prevalent") %>%
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

scores_sum_pheno <- scores_sum %>%
   left_join(all_ids, by = c("id"="eid")) %>%
   filter(!is.na(pheno)) %>% 
   left_join(fam_hist, by = c("id"="eid")) %>% 
   mutate(ill = ifelse(illness == "Breast cancer", 1, 0 )) %>% 
   mutate(ill = replace_na(ill,0)) %>% 
   left_join(ukbb_pc_df %>% dplyr::select(eid, genetic_principal_components_f22009_0_1,
                                           genetic_principal_components_f22009_0_2  ),
              by = c("id" ="eid"))

scores_sum_pheno %>% 
   group_by(pheno) %>% 
   count(illness)
#scores_sum_pheno %>% 
#  mutate(FID = id, IID = id) %>% 
#   select(FID,IID,pheno) -> target_pheno
#fwrite(target_pheno, file = file.path(export_path,"ukbb_pheno_test_output","target_pheno.txt"), 
#       sep="\t", row.names = FALSE)
#scores_sum_pheno %>% 
#   mutate(FID = id, IID = id) %>% 
#   select(FID,IID,age_when_attended_assessment_centre_f21003_0_0 ,
#          genetic_principal_components_f22009_0_1,
#          genetic_principal_components_f22009_0_2) -> cov

#fwrite( cov, file = file.path(export_path,"ukbb_pheno_test_output","cov.txt"), 
#       sep="\t", row.names = FALSE)
# 4.0 Exploratory plots ----

hist(scores_sum_pheno$sum)

score_interest %>%
   ggplot(aes(x=scale(exome_score_sum),after_stat(density))) + 
   geom_histogram(alpha=0.3) +
   geom_density(alpha=0.7) +
   theme_pubr()



score_interest %>%
   ggplot(aes(x=as.factor(pheno), y=exome_score_sum, fill=as.factor(pheno))) +
   geom_boxplot() +
   theme_pubr()

t.test( x= scores_sum_pheno %>% filter(pheno == 1) %>% dplyr::select(sum),
        y= scores_sum_pheno %>% filter(pheno == 0)%>% dplyr::select(sum))

scores_sum %>%
   left_join(ukbb_pc_df %>% dplyr::select(eid, genetic_principal_components_f22009_0_1,
   genetic_principal_components_f22009_0_2,
   genetic_ethnic_grouping_f22006_0_0),
   by = c("id" ="eid")) %>%
   mutate(eth = case_when(
   genetic_ethnic_grouping_f22006_0_0 == "Caucasian" ~ 1,
   TRUE ~ 0
   )) %>%
   dplyr::select(sum,eth) -> eth_score

eth_score %>%
   ggplot(aes(x=as.factor(eth), y=sum, fill=as.factor(eth))) +
   geom_boxplot() +
# geom_point(aes(colour = pheno)) +
   theme_pubr()

t.test( x= eth_score %>% filter(eth == 1) %>% select(sum),
        y= eth_score %>% filter(eth == 0)%>% select(sum))

ukbb_pc_df %>%
   ggplot(aes(genetic_principal_components_f22009_0_1, genetic_principal_components_f22009_0_2,
         color= ethnic_background_f21000_0_0)) +
   geom_point()+
   theme_pubr()



# 5.0 Performing logistic regression ----

scores_sum_pheno %>%
   left_join(ukbb_pc_df %>% dplyr::select(eid, genetic_principal_components_f22009_0_1,
                                   genetic_principal_components_f22009_0_2,
                                   age_when_attended_assessment_centre_f21003_0_0  ),
             by = c("id" ="eid")) -> score_pc
# association model


model <- glm(pheno ~  
                sum + age + ill +
                
                genetic_principal_components_f22009_0_1 +
                genetic_principal_components_f22009_0_2, 
             scores_sum_pheno , family = "binomial") 

tidy(model)
summary(model) -> sum_plink_log

sum_plink_log$coefficients %>% as.data.frame() %>% str

# calculating Nagelkerke R2

nagelkerke(model)
lrm(model)

nagelkerke_r2 <- nagelkerke(model)$Pseudo.R.squared.for.model.vs.null[3]

# calculating odd ratio
logistic.display(model, decimal = 3,crude = FALSE,simplified = FALSE)

# training model

inTrain<- createDataPartition(y=scores_sum_pheno$pheno, p=0.75, list=FALSE)
training<-scores_sum_pheno[inTrain,]

testing<-scores_sum_pheno[-inTrain,]

glm(pheno ~  
       sum + 
       age +
       ill +
       
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2, 
    training , family = "binomial") -> model



predict(model, testing) -> prob

testing$prob=prob

# calculating AUC


g <- roc(pheno ~ prob, data = testing, plot=TRUE,
         print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

library(plotROC)
rocplot <- ggplot(testing, aes(m = prob, d = pheno ))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot +
   style_roc(guide = FALSE) +
   geom_rocci()+
   annotate("text", x = .75, y = .25,
            label = paste("AUC =", round(calc_auc(rocplot)$AUC, 2))) +
   theme_pubr()


#quantile(scores_sum_pheno$sum, prob = seq(0, 1, length = 11), type = 5)

# prevalence plots ----
score_interest %>% 
   mutate(percentile = ntile(score_interest$exome_score_sum,100) ) %>%
      ggplot(aes(x=as.factor(pheno), y=percentile, fill=as.factor(pheno))) +
      geom_boxplot() +
      scale_y_continuous(breaks = seq(0, 100, by = 10))+
      theme_pubr(legend = "none")

score_interest %>% 
   mutate(percentile = ntile(score_interest$exome_score_sum,5) ) %>% 
   mutate(pheno = ifelse(pheno == "yes" , 1, 0)) %>% 
   group_by(percentile) %>% 
   summarise(prevalence = 100*sum(pheno)/n()) %>% 
   ggplot(aes(x= percentile, y=prevalence,color = percentile))+
   scale_color_gradient(low = "#56B1F7", high = "#132B43")+
   geom_point()+
   theme_pubr()

# Odd ratio ----

cut_list <- seq.int(from = median(scores_sum_pheno$sum) , to = round(max(scores_sum_pheno$sum),3), 
                    length.out = 100 )


odd_ratio <- function(data, threshold){
   final <- data.frame(matrix(ncol=4, nrow=length(threshold)))
   x <- c("cut_thre", "OR", "CI","prev_reminder")
   colnames(final) <- x
   for (i in 1:length(threshold))
      {
   data %>%  
   mutate(thre = ifelse(sum > threshold[i] ,"above","below")) %>%
   group_by(pheno, thre) %>% 
   summarise(total= n()) %>%
   spread(pheno,total) %>%
   mutate(prev_reminder = round(100*(`1`+`0`)/nrow(data),3)) -> tmp 
   
   tmp %>%       
   column_to_rownames("thre") %>% 
   dplyr::select(`1`,`0`) %>% 
   fisher.test() -> e

final[i,] <- 
   tribble(
      ~cut_thre,   ~OR,  ~CI,  ~prev_reminder,
       threshold[i],  as.numeric(e$estimate), toString(e$conf.int), paste(tmp$prev_reminder[1],"%") )
   }

final %>% 
   mutate(floor_OR = floor(OR)) %>% 
   group_by(floor_OR) %>% 
   slice(which.min(OR))
}

odd_ratio(scores_sum_pheno,cut_list)


dat <- with(density(scores_sum_pheno$sum), data.frame(x, y))
ggplot(data = dat, mapping = aes(x = x, y = y)) +
   geom_line()+
   geom_vline(xintercept=0.0303, size=0.5, color="#69baf8")+
   geom_area(mapping = aes(x = ifelse(x>0.0303 & x< 0.0787 , x, 0)), fill = "#69baf8") +
   geom_vline(xintercept=0.0787, size=0.5, color="#0862a8")+
   geom_area(mapping = aes(x = ifelse(x>=0.0787 & x< 0.0982 , x, 0)), fill = "#0862a8") +
   geom_vline(xintercept=0.0982, size=0.5, color="#132B43")+
   geom_area(mapping = aes(x = ifelse(x>=0.0982 , x, 0)), fill = "#0862a8") +
   ylim(0, 18) +
   theme_pubr()


# Odd Ratio scale ----
cut_list <- seq.int(from = 0.01 , to = 4.7, length.out = 100 )

odd_ratio <- function(data, threshold){
   final <- data.frame(matrix(ncol=4, nrow=length(threshold)))
   x <- c("cut_thre", "OR", "CI","prev_reminder")
   colnames(final) <- x
   for (i in 1:length(threshold))
   {
      data %>%  
         mutate(thre = ifelse(scale(sum) > threshold[i] ,"above","below")) %>%
         group_by(pheno, thre) %>% 
         summarise(total= n()) %>%
         spread(pheno,total) %>%
         mutate(prev_reminder = round(100*(`1`+`0`)/nrow(scores_sum_pheno),3)) -> tmp 
      
      tmp %>%       
         column_to_rownames("thre") %>% 
         dplyr::select(`1`,`0`) %>% 
         fisher.test() -> e
      
      final[i,] <- 
         tribble(
            ~cut_thre,   ~OR,  ~CI,  ~prev_reminder,
            threshold[i],  as.numeric(e$estimate), toString(e$conf.int), paste(tmp$prev_reminder[1],"%") )
   }
   
   final %>% 
      mutate(floor_OR = floor(OR)) %>% 
      group_by(floor_OR) %>% 
      slice(which.min(OR))
}

odd_ratio(scores_sum_pheno,cut_list)
   
dat <- with(density(scale(scores_sum_pheno$sum)), data.frame(x, y))
ggplot(data = dat, mapping = aes(x = x, y = y)) +
   geom_line()+
   geom_vline(xintercept=0.721, size=0.5, color="#69baf8")+
   geom_area(mapping = aes(x = ifelse(x>0.721 & x< 2.71 , x, 0)), fill = "#69baf8") +
   geom_vline(xintercept=2.71, size=0.5, color="#0862a8")+
   geom_area(mapping = aes(x = ifelse(x>=2.71 & x< 3.52 , x, 0)), fill = "#0862a8") +
   geom_vline(xintercept=3.52, size=0.5, color="#132B43")+
   geom_area(mapping = aes(x = ifelse(x>=3.52 , x, 0)), fill = "#0862a8") +
   ylim(0, 0.4) +
   theme_pubr()

scores_sum_pheno %>%  
   mutate(thre = ifelse(scale(sum) > 4.7000000 ,"above","below")) %>%
   group_by(pheno, thre) %>% 
   summarise(total= n()) %>%
   spread(pheno,total) %>%
   mutate(prev_reminder = round(100*(`1`+`0`)/nrow(scores_sum_pheno),3)) -> tmp 

tmp %>%       
   column_to_rownames("thre") %>% 
   dplyr::select(`1`,`0`) %>% 
   fisher.test() 




scores_sum_pheno %>%  
   mutate(thre = ifelse(sum > 0.05 ,"above","below")) %>%
   group_by(pheno, thre) %>% 
   summarise(total= n()) %>%
   spread(pheno,total) %>%
   mutate(prev_reminder = round(100*(`1`+`0`)/nrow(scores_sum_pheno),3)) -> tmp 

tmp %>%       
   column_to_rownames("thre") %>% 
   dplyr::select(`1`,`0`) %>% 
   fisher.test() 


fread("/home/opc/storage/313_snps/gwas-association-downloaded_2020-07-16-pubmedId_29059683.tsv") ->eee

# PRScise results ----
## import files 
prsice_0.2_best <- fread("/home/opc/storage/output/output/score_bcac/breast_fem_british_bcac_0.2.best") 

scores_sum_pheno %>% 
   inner_join(prsice_0.2_best, by = c("id"="IID")) %>% 
   ggplot(aes(scale(sum),scale(PRS))) +
   geom_point()

scores_sum_pheno %>%
   inner_join(prsice_0.2_best, by = c("id"="IID")) %>% 
   ggplot(aes(x=as.factor(pheno), y=PRS, fill=as.factor(pheno))) +
   geom_boxplot() +
   theme_pubr()

scores_sum_pheno %>% 
   inner_join(prsice_0.2_best, by = c("id"="IID")) -> mix_prs

cor.test(mix_prs$sum,mix_prs$PRS)


model <- glm(pheno ~ PRS + age + ill +

                
                genetic_principal_components_f22009_0_1 +
                genetic_principal_components_f22009_0_2, 
             mix_prs , family = "binomial") 

summary(model) 
tidy(model)

inTrain<- createDataPartition(y=mix_prs$pheno, p=0.75, list=FALSE)
training<-mix_prs[inTrain,]

testing<-mix_prs[-inTrain,]

glm(pheno ~  
       sum + age + ill +
       
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2, 
    training , family = "binomial") -> model

glm(pheno ~ sum + genetic_principal_components_f22009_0_1 + genetic_principal_components_f22009_0_2 ,
    data=training,
    
    family = "binomial") -> model

predict(model, testing) -> prob

testing$prob=prob

# calculating AUC


g <- roc(pheno ~ prob, data = testing, plot=TRUE,
         print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

nagelkerke_r2 <- nagelkerke(model)$Pseudo.R.squared.for.model.vs.null[3]

scores_sum_pheno %>%
   inner_join(prsice_0.2_best, by = c("id"="IID")) -> sum_prs

# Lassossum ----
install.packages("devtools")
install_github("tshmak/lassosum")
library(lassosum)


# import files

## summary statistics
sum_stats <- fread("/home/opc/storage/output/output/bcac_sumstats/bcac_varid_matched_summary_stats.txt") %>% 
   filter(bcac_gwas_all_beta != "NULL" & bcac_gwas_all_P1df != "NULL") 

### derive SNP-wise correlations.
cor <- p2cor(p = as.numeric(sum_stats$bcac_gwas_all_P1df), n = nrow(sum_stats),
             sign=as.numeric(sum_stats$bcac_gwas_all_beta)) #n=sample size sumstats

## phenotype and covariates

target_pheno <- fread("/home/opc/storage/output/ukbb_pheno_test_output/target_pheno.txt")

cov <-  fread("/home/opc/storage/output/ukbb_pheno_test_output/cov.txt")


#This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
# Other alternatives available. Type ?lassosum.pipeline for more details. 
LDblocks <- "EUR.hg19"

# open a cluster for parallel computation and define and output list
cluster=makeCluster(10, type="FORK")
#output=list()
no_cores <- availableCores() - 1
cluster2 <- makeClusterPSOCK(workers = 10)
# run lassosum by chromosome

ref_path <- "/home/opc/storage/output/Lassosum/CEU_1KG/"
test_path <- "/home/opc/storage/uk_biobank/genotype_data/array/"

options(future.globals.maxSize= 2000000000)

no_cores <- availableCores() - 1
plan(multisession, workers = no_cores)

laso_out <- c(1:22) %>% 
   future_map_dfr(function(x) {
      lassosum.pipeline(cor=cor, chr=sum_stats$CHR, pos=sum_stats$BP, 
                        A1=sum_stats$A0, A2=sum_stats$A1,
                        ref.bfile=paste(ref_path,x,"_CEU1kG_ref",sep=""), 
                        test.bfile=paste(test_path, "ukb_chr",x,sep=""),
                        LDblocks = LDblocks)}) 


output =lapply(c(1:22),
              function(x) {lassosum.pipeline(cor=cor, chr=sum_stats$CHR, pos=sum_stats$BP, 
                                             A1=sum_stats$A0, A2=sum_stats$A1,
                                             ref.bfile=paste(ref_path,x,"_CEU1kG_ref",sep=""), 
                                             test.bfile=paste(test_path, "ukb_chr",x,sep=""),
                                             LDblocks = LDblocks,cluster = cluster)}) 

lassosum.pipeline(cor=cor, chr=sum_stats$CHR, pos=sum_stats$BP, 
                  A1=sum_stats$A0, A2=sum_stats$A1,
                  ref.bfile=paste(ref_path,1,"_CEU1kG_ref",sep=""), 
                  test.bfile=paste(test_path, "ukb_chr",1,sep=""),
                  LDblocks = LDblocks)
