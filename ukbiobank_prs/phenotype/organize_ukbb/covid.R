# 1.0 Load libraries ----
print("install and loading dependencies packages")


installed.packages("pacman") 

library("pacman")

#devtools::install_github("kenhanscombe/ukbtools", build_vignettes = TRUE, dependencies = TRUE, force = TRUE)

pacman::p_load(tidyverse     # tidy daya
               ,
               ukbtools    # workhorse package for organizing ukbb
               ,
               tictoc      # checking time
               ,
               furrr       # functional programing with map()
               ,
               future      # parallel computing
               ,
               arrow       # apache arrow provides powerful package for sharing data
               ,
               data.table  # exporting data as .txt
               ,
               ggpubr
               , 
               lubridate,
               eeptools)     # add themes and visualization prop. to ggplot   

# 2.0 Importing Files ----

print("source change_me.r file")
source(here::here("phenotype", "organize_ukbb", "change_me.R"))

fread(file = file.path("/home/opc/storage/uk_biobank/phenotype_data/ukbiobank_covid" ,
                       'covid19_result.txt')) -> covid_ukb

fread(file = file.path("/home/opc/storage/uk_biobank/phenotype_data/ukbiobank_covid" ,
                       'death.txt')) -> death
fread(file = file.path("/home/opc/storage/uk_biobank/phenotype_data/ukbiobank_covid" ,
                       'death_cause.txt')) -> death_cause

covid_ukb %>% 
   left_join(death) %>% 
   mutate(severity = ifelse( result == 1 & origin == 1 , 1 , 0 ),
          death = ifelse(is.na(date_of_death) == FALSE, 1,0 )) -> covid_ukb 

covid_ukb %>% select(eid,result,severity,death) %>% 
   left_join(ukbb_df %>% select(eid,sex_f31_0_0,matches("genetic_principal_components")), by = c("eid"="eid") ) %>%
   left_join(death %>% select(eid,date_of_death)) %>% 
   left_join(birth_ukbb_df) %>% 
   mutate(current_date = as_date(now()),
          date_of_death = as.Date(date_of_death, "%m/%d/%y"),
          current_date = as.Date(ifelse(death == 1 , date_of_death ,current_date), origin = "1970-01-01"),
          current_age = time_length(difftime(current_date, date_of_birth), "years")) %>% 
   select(-date_of_death,-date_of_birth,-current_date) -> covid_ukb_final

fwrite(covid_ukb_final , file = file.path(export_path, "output","ukbb_pheno_test_output",
                                 "covid_ukbb.txt"), 
       sep="\t", row.names = FALSE)







covid_ukb  %>% 
   left_join(death_cause, by = "eid") %>% 
   filter(severity == 1 & death == 1) %>% count(cause_icd10) %>% 
   arrange(desc(n)) %>% 
   left_join(icd10codes, c("cause_icd10" = "code"))
   
   



