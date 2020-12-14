dir_exome_scores <- "/home/opc/storage/output/output/Exome_scores/200k"
no_cores <- availableCores() - 1

plan(multisession, workers = no_cores)

tic()
genepy_exome_scores_1 <- list.files(dir_exome_scores , full.names = T) %>%
   future_map(fread) %>% 
   map(~ subset(.x,select=which(!duplicated(names(.x))) )) %>% 
   reduce(full_join, by= c("IID" = "IID"))
toc()


tic()
feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
                          "all_chr_genepy_score.feather")
arrow::write_feather(genepy_exome_scores_1 , feather_file)
toc()

tic()
feather_file <- file.path("/home/opc/storage/output/output/Exome_scores",
                          "prs_gene_score_cad.feather")
prs_gene_score_cad <- arrow::read_feather(feather_file) %>% as.data.frame()
toc()
