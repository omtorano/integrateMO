## code to prepare `test_RRBS_feature_map` dataset goes here
## Data import
test_RRBS_feature_map <- read.csv("C:/Users/otorano/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/EE2/test_rrbs_package/ML-0_1000-1000_0.8_Unique_Features_to_Genes.csv")
saveRDS(test_RRBS_feature_map, "data-raw/test_RRBS_feature_map.RDS")
usethis::use_data(test_RRBS_feature_map, overwrite = TRUE)
