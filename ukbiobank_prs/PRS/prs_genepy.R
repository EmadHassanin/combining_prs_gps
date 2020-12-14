exome_scores <- fread("/home/opc/storage/output/output/Exome_scores/Genome_Wide_Exome_Score")

exome_scores %>% colnames() %>% make.unique() -> new_col_names
colnames(exome_scores) <- new_col_names
#exome_scores %>% dplyr::select(IID,BRCA1,BRCA2,PALB2,CHEK2,PTEN,CDH1) -> exome_scores_breast
#exome_scores_sum <- exome_scores %>% 
#   mutate( exome_score_sum =  rowSums(exome_scores %>% dplyr::select(- IID), na.rm = TRUE))

conversion_table <- fread("/home/opc/storage/output/output/Exome_scores/conversion_tables")


exome_scores %>%  
   inner_join(conversion_table %>% dplyr::select(V1,V7), by = c('IID' = 'V1')) %>% 
   inner_join(prev_scores_sum_313, by = c('V7' = 'id')) -> prs_gene_score

exome_scores %>%  
   inner_join(conversion_table %>% dplyr::select(V1,V7), by = c('IID' = 'V1')) %>% 
   inner_join(pheno_icd_registered, by = c('V7' = 'eid')) %>% 
   group_by(ethnic_background_f21000_0_0,pheno) %>% 
   count() %>% 
   view()


options(expressions = 5e5)


out_term <- data_frame(matrix(ncol=2, nrow=length(seq(from=3, to=30, by=2212))))

model2 <- map(seq(from=1, to=17608, by= 1),
             possibly(function(x) {
         glm(pheno ~ . - IID - V7   , 
             
             prs_gene_score %>% dplyr::select(x, IID, pheno, V7, genetic_principal_components_f22009_0_1,
                                              genetic_principal_components_f22009_0_1, age), 
             family = "binomial") %>% 
             tidy()}, NA_real_) )  %>% 
   reduce(rbind) %>%
   distinct(term, .keep_all = TRUE)

fwrite(model2, file = file.path("/home/opc/storage/output/output/" ,
                                'model.txt'), 
           sep="\t", row.names = FALSE)
fread(file = file.path("/home/opc/storage/output/output/Exome_scores" ,
                       'model.txt')) -> model5

model <- model5 %>% 
   filter(str_starts(term, "\\w+" ), !is.na(p.value), p.value < 0.05)

prs_gene_score %>% 
   select(model$term,pheno,IID,genetic_principal_components_f22009_0_1, 
          genetic_principal_components_f22009_0_2)  -> score_interest
prs_gene_score %>% 
mutate( exome_score_sum =  rowSums(prs_gene_score %>% dplyr::select(-IID, -age,-sum, -ill, -illness,
                        -pheno,-genetic_principal_components_f22009_0_1,
                       -genetic_principal_components_f22009_0_2), na.rm = TRUE)) ->score_interest2

glm(pheno ~ BRCA1+BRCA2+PALB2+CHEK2+PTEN+CDH1 + sum + age +
       #ill +
       
       genetic_principal_components_f22009_0_1 +
       genetic_principal_components_f22009_0_2,
    prs_gene_score,
       
    family = "binomial") -> model3
prs_gene_score %>% dplyr::select(one_of(genes) , sum, pheno, age, genetic_principal_components_f22009_0_1, 
                                 genetic_principal_components_f22009_0_2)  -> prs_gene_score_2
summary(model3)
set.seed(1)
inTrain<- createDataPartition(y=score_interest2$pheno, p=0.75, times = 10,list=FALSE)
training<-score_interest2[inTrain,]

testing<-score_interest2[-inTrain,]

glm(pheno ~  exome_score_sum + sum + age,

    data=training, family = "binomial") -> model_1

glm(pheno ~ sum + age  ,
    #      genetic_principal_components_f22009_0_1 ,
    #     genetic_principal_components_f22009_0_2,
    data=training, family = "binomial") -> model_2

glm(pheno ~  age ,
    #ill +
    
    #  genetic_principal_components_f22009_0_1 +
    #    genetic_principal_components_f22009_0_2,
    data=training, family = "binomial") -> model_3



predict(model_1, testing) -> prob_1
predict(model_2, testing) -> prob_2
predict(model_3, testing) -> prob_3
testing$prob_1=prob_1
testing$prob_2=prob_2
testing$prob_3=prob_3
# calculating AUC


g_1 <- roc(pheno ~ prob_1, data = testing, plot=TRUE,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

g_2 <- roc(pheno ~ prob_2, data = testing, plot=TRUE,
         print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

g_3 <- roc(pheno ~ prob_3, data = testing, plot=TRUE,
           print.auc=TRUE,ci=TRUE, legacy.axes = TRUE)

rocplot <- ggplot(testing, aes(m = prob_1, d = pheno ))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot <- rocplot +
   style_roc(guide = FALSE) +
   geom_rocci()+
   annotate("text", x = .75, y = .25,
            label = paste("AUC =", round(calc_auc(rocplot)$AUC, 2))) +
   theme_pubr()

g <- ggplot(testing, aes(d = pheno , color = name))
g <- g + geom_roc(aes(m=prob_1), colour="red",n.cuts=20,labels=FALSE) 
g <- g+ geom_roc(aes(m=prob_2), colour="green",n.cuts=20,labels=FALSE)
g <- g+ geom_roc(aes(m=prob_3), colour="blue",n.cuts=20,labels=FALSE) 
g +
   annotate("text", x = .75, y = .25,
            label = paste("AUC =", round(calc_auc(g)$AUC, 2))) +
   theme_pubr()

longtest <- melt_roc(testing, "pheno", c("prob_1", "prob_2","prob_3")) %>% 
   mutate(name = case_when(
      name == "prob_1" ~ "genepy + prs", 
      name == "prob_2" ~"prs",
      name == "prob_3" ~"genepy"))
ggplot(longtest, aes(d = D, m = M, color = name)) + 
   geom_roc(n.cuts=20,labels=FALSE) ->p
direct_label(p,labels = c("prs + genepy" ,"prs")) + theme_pubr()

p +
   annotation_custom(tableGrob(aucs, rows=NULL, theme=ttheme_minimal()), 
                     xmin=0.5, xmax=1, ymin=0, ymax=0.5) + 
   theme_pubr()

calc_auc(p)

model <- unique(longtest$name)
model_info <- data.frame(model,
                         group = rank(model))
left_join(model_info, calc_auc(p)) %>%
   select(-group, -PANEL) %>%
   arrange(desc(AUC)) %>% 
   mutate( AUC = round(AUC,3)) -> aucs


genes <- read_excel("/home/opc/storage/output/output/Exome_scores/41467_2018_3411_MOESM10_ESM.xlsx", 
                    skip = 1) %>% select(1)
genes$`RefSeq gene` -> genes
genes %>% append(c("BRCA1","BRCA2","PALB2","CHEK2","PTEN","CDH1")) -> genes
c("BRCA1","BRCA2","PALB2","CHEK2","PTEN","CDH1") -> genes
genes <- c("AR","ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CASP8", "CDH1", "CHEK2", "CTLA4","COMT" ,"CYP19A1","CYP1A1",
            "CYP1B1","CYP2D6","CYP17","CYP19","ER","FGFR2", "H19","APC", "HER2",
  "LSP1", "MAP3K1", "MRE11", "NBN", "PALB2", "PTEN","PR" ,"RAD51", "RAD51C", "STK11", "TERT", "TOX3", "TP53", 
  "XRCC2", "XRCC3", "LKB1" , "HRAS1" , "NAT1" , "NAT2" , "GSTM1", "GSTP1" , "GSTT1", "UGT1A1" ,"APOE", "CYP2E1")
prs_gene_score %>% 
   dplyr::select(one_of(genes),sum , pheno, age,IID,
                 genetic_principal_components_f22009_0_1, 
                 genetic_principal_components_f22009_0_2)  -> prs_gene_score_2
summary(model3)



