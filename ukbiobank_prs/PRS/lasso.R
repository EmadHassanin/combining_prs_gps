
pacman::p_load(tidyverse,data.table, tictoc, furrr, future,
               lassosum,
               parallel) 

#install_github("tshmak/lassosum@v0.4.4")

# import files

## summary statistics
write("TMP = /home/opc/storage/output/output/tmp", file=file.path('~/.Renviron'))

sum_stats <- fread("/home/opc/storage/output/output/bcac_sumstats/bcac_varid_matched_summary_stats.txt") %>% 
   filter(bcac_gwas_all_beta != "NULL" & bcac_gwas_all_P1df != "NULL") 
sum_stats %>% 
   mutate(bcac_gwas_all_P1d = as.numeric(bcac_gwas_all_P1df),
          bcac_gwas_all_beta = as.numeric(bcac_gwas_all_beta)
          ) -> sum_stats

#fwrite(sum_stats, file = "/home/opc/storage/output/output/bcac_sumstats/bcac_sum.txt", 
#           sep="\t", row.names = FALSE)

### derive SNP-wise correlations.
cor <- p2cor(p = as.numeric(sum_stats$bcac_gwas_all_P1df), n = nrow(sum_stats),
             sign=as.numeric(sum_stats$bcac_gwas_all_beta)) #n=sample size sumstats

## phenotype and covariates

#target_pheno <- fread("/home/opc/storage/output/ukbb_pheno_test_output/target_pheno.txt")

#cov <-  fread("/home/opc/storage/output/ukbb_pheno_test_output/cov.txt")


#This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
# Other alternatives available. Type ?lassosum.pipeline for more details. 
LDblocks <- "EUR.hg19"

# open a cluster for parallel computation and define and output list
cluster=makeCluster(40, type="FORK")
#output=list()
#no_cores <- availableCores() - 1
#cluster2 <- makeClusterPSOCK(workers = 10)
# run lassosum by chromosome

#ref_path <- "/home/opc/storage/output/Lassosum/CEU_1KG/"
test_path <- "/home/opc/storage/uk_biobank/genotype_data/array/"

options(future.globals.maxSize= 2000000000)

#no_cores <- availableCores() - 1
#plan(multisession, workers = no_cores)

# tic()
# laso_out <-
#    lapply( c(1:22),
#            function(x) {
#       lassosum.pipeline(cor=cor, chr=sum_stats$CHR, pos=sum_stats$BP, 
#                         A1=sum_stats$A0, A2=sum_stats$A1,
#                         test.bfile=paste(test_path, "ukb_chr",x,sep=""),
#                         LDblocks = LDblocks,cluster = cluster)}) 
# toc()
# save(laso_out, file=file.path(export_path,"lasso.RData"))

#output =sapply(c(1:22),
#               function(x) {lassosum.pipeline(cor=cor, chr=sum_stats$CHR, pos=sum_stats$BP, 
#                                              A1=sum_stats$A0, A2=sum_stats$A1,
#                                              ref.bfile=paste(ref_path,x,"_CEU1kG_ref",sep=""), 
#                                              test.bfile=paste(test_path, "ukb_chr",x,sep=""),
#                                              LDblocks = LDblocks,cluster = cluster)}) 

tic()
laso_out <-lassosum.pipeline(cor=cor, chr=sum_stats$CHR, pos=sum_stats$BP, snp=sum_stats$SNP, 
                                A1=sum_stats$A0, A2=sum_stats$A1,sample = 5000,
                                test.bfile=file.path(test_path,"ukb_chr22"),
                                LDblocks = LDblocks,cluster = cluster) 
toc()

exp_path <- "/home/opc/storage/output/output/lasso_score"
save(laso_out, file=file.path(exp_path ,"lasso2.RData"))
