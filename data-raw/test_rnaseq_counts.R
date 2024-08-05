## code to prepare `test_RNAseq_counts` dataset goes here
## Data import
# rrbs data
liver_rrbs <- readRDS("C:/Users/otorano/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/EE2/test_rrbs_package/ML-0_1000-1000_0.8_Mval.RDS")
colnames(liver_rrbs) <- gsub("_Mval","",colnames(liver_rrbs))
# gene level counts data
all_counts <- read.delim("C:/Users/otorano/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/EE2/RNAseq/rsem.merged.gene_counts.tsv")
rownames(all_counts) <- all_counts$gene_id
all_counts <- all_counts[, -c(1, 2)]
colnames(all_counts) <- gsub("\\.", "-", colnames(all_counts))
colnames(all_counts) <- gsub("2-5", "2.5", colnames(all_counts))
# metadata
ee2_meta <- read.delim("C:/Users/OTORANO/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/EE2/RNAseq/ee2_meta.txt")
ee2_meta <- dplyr::distinct(ee2_meta)
rownames(ee2_meta) <- ee2_meta$sample
match_samp_ee2 <- table(c(colnames(all_counts), colnames(liver_rrbs), ee2_meta$sample))
match_samp_ee2 <- names(match_samp_ee2[match_samp_ee2 > 2])
liver_meta <- ee2_meta[ee2_meta$sample%in%match_samp_ee2, ]
colnames(liver_meta)[2] <- "TRT"
liver_counts <- all_counts[, liver_meta$sample]
# reorder rrbs rows to match meta and counts data
liver_rrbs <- liver_rrbs[, match(liver_meta$sample, colnames(liver_rrbs))]
# subset top 10000 features associated with WGCNA modules that separate counts by treatment & (+) correlated methylation modules
test_rnaseq_counts <- liver_counts
genes_to_include <- read.csv("C:/Users/otorano/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/EE2/integrateMO_WGCNA_2023-11-03 15.10.28/Module_membership_rnaseq_counts.csv",
                             row.names = 1)
keep <- rownames(genes_to_include[order(genes_to_include$MEcyan, genes_to_include$MEgreenyellow,
                                        decreasing = TRUE), ])[1:500]
test_rnaseq_counts <- liver_counts[keep, ]
saveRDS(test_rnaseq_counts, "data-raw/test_rnaseq_counts.RDS")
usethis::use_data(test_rnaseq_counts, overwrite = TRUE)
