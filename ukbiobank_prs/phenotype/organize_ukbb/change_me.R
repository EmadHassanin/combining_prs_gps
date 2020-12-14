# Config for potential phenotypes data from ukbiobank ----

# 1.0 raw data directory ----
#' this directory should contain .tab, .r, and .html files

data_path <- "/home/opc/storage/uk_biobank/phenotype_data/ukb_df"
file_name <- "ukb_df.rds"
icd_file_name <- "icdcodes_all.rds"
fam_files_path <- "/home/opc/storage/output/fam_files"

# 2.0 identify icd codes ----
# identify both icd versions 9 and 10. 
# if only one version needed, leave the other empty (c(""))

identify_icd <- 
tribble(
  ~phenotype      , ~icd9         , ~icd10,
  "breast_cancer" , "^(174[0-9])" , "^(C50[0-9])",
  "cad"           , "^(410[0-9]|411[0-9]|412[0-9]|42979)","^(I21[0-9]|I22[0-9]|I23[0-9]|I241|I252)",
  "diabetes_2"    , "no"            , "^(E11[0-9])"
)

        


col_vars <- 
tribble(
  ~col_var             , ~col_name,
  "id"                 , "eid",
  "sex"                , "sex_f31_0_0",
  "age_at_assessment"  , "age_when_attended_assessment_centre_f21003_0_0",
  "age_at_death"       , "age_at_death_f40007_2_0",
  "body_mass_index"    , "body_mass_index_bmi_f21001_0_0",
  "socioeconomic"      , "townsend_deprivation_index_at_recruitment_f189_0_0",
  "ethnicity"          , "ethnic_background_f21000_0_0",
  "employment"         , "current_employment_status_f6142_0_0",
  "centre"             , "uk_biobank_assessment_centre_f54_0_0",
  "Genotype_batch"     , "genotype_measurement_batch_f22000_0_0",
  "Genotype_plate"     , "genotype_measurement_plate_f22007_0_0",
  "Genotype_well"      , "genotype_measurement_well_f22008_0_0" ,
  "Genetic_sex"        , "genetic_sex_f22001_0_0",
  "Genetic_kinship"    , "genetic_kinship_to_other_participants_f22021_0_0",
  "Genetic_ethnic"     , "genetic_ethnic_grouping_f22006_0_0",
  "Sex chromosome aneuploidy", "sex_chromosome_aneuploidy_f22019_0_0",
  "Outliers for heterozygosity or missing rate", "outliers_for_heterozygosity_or_missing_rate_f22027_0_0",
  "Heterozygosity"     , "heterozygosity_f22003_0_0",
  "Heterozygosity, PCA corrected", "heterozygosity_pca_corrected_f22004_0_0",
  "Missingness"        , "missingness_f22005_0_0" ,
  #"Genetic principal components", "" ,
  "Used in genetic principal components", "used_in_genetic_principal_components_f22020_0_0"
)

# 4.0 export data directory ----


export_path <- "/home/opc/storage/output"





