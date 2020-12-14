# 1.0 Load libraries ----

library("pacman")

#devtools::install_github("kenhanscombe/ukbtools", build_vignettes = TRUE, dependencies = TRUE, force = TRUE)


pacman::p_load(tidyverse, ukbtools, tictoc, furrr, future, data.table)
# 2.0 import data ----
dir_path <- "/home/opc/storage/output/snp_list/khera_breast_snps"
dir_ukbiobank <- "/home/opc/storage/uk_biobank/genotype_data/array/"
SNP_khera <- fread(file.path(dir_path,"BreastCancer_PRS_PT_r2_0.2_p_0.0005_v3.txt")) %>% 
   mutate(chr_pos = paste0(chr ,":",position_hg19)) 
   
#SNP_matched <- fread(paste0(dir_ukbiobank,"ukb_chr3.bim"))

#SNP_matched <- SNP_matched %>% 
#   mutate(chr_pos = paste0(V1,":",V4))

no_cores <- availableCores() - 1
plan(multisession, workers = no_cores)

tic()

c(1:22) %>% 
   future_map_dfr(function(x) {
SNP_khera %>% 
   inner_join(
      fread(paste0(dir_ukbiobank,"ukb_chr",x,".bim")) %>% 
         mutate(chr_pos = paste0(V1,":",V4)) , by= c("chr_pos" = "chr_pos")
   ) } ) -> snps_khera_final

snps_khera_final %>% 
   dplyr::select(V2,effect_allele, effect_weight) %>% 
   unique() -> snps_khera_final_2

fwrite(snps_khera_final, file = file.path(dir_path,"snps_khera_final.txt"), 
           sep="\t", row.names = FALSE)


#colnames(snps_313_final) <- NULL
#colnames(snps_313_final) <- c("snp_id","effect_allele","overall_breast_cancer")

fread("/home/opc/storage/output/bcac_sumstats/bcac_varid_matched_summary_stats.txt") -> bcac_sum

library(qqman)
library("QCEWAS")

bcac_sum %>% 
   filter(bcac_gwas_all_P1df < 0.005) %>% 
   mutate(bcac_gwas_all_P1df = as.numeric(bcac_gwas_all_P1df))-> bcac_filtered

   
manhattan(bcac_filtered,  p = "bcac_gwas_all_P1df")
sample_n(bcac_sum, 100000) %>% dplyr::select( bcac_gwas_all_P1df) %>%  as.vector()-> bcac_pvalue
qq(as.numeric(bcac_pvalue$bcac_gwas_all_P1df))

P_lambda(as.numeric(bcac_pvalue$bcac_gwas_all_P1df))
hist(as.numeric(bcac_filtered$bcac_gwas_all_beta))

plot(abs(as.numeric(bcac_filtered$bcac_gwas_all_beta)), -log10(as.numeric(bcac_filtered$bcac_gwas_all_P1df)))


bcac_filtered %>% filter(bcac_gwas_all_P1df < 0.00000001) %>% select(bcac_gwas_all_beta ) -> beta_filtered


hist(as.numeric(beta_filtered$bcac_gwas_all_beta), 20)

