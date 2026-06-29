library(data.table)
library(dplyr)
library(stringr)
library(argparse)

#set up argparse
parser <- ArgumentParser()
parser$add_argument("--cell_type", help="Single-Cell Type")
args <- parser$parse_args()

#build file read in patterns - updated for new naming convention
pattern <- paste0("~/sc_pwas_in_ukb/output/", args$cell_type, "_cell/*.csv")
files <- Sys.glob(pattern)
cat("Number of files found:", length(files), "\n\n")

#initialize empty list
all_results <- list()
all_results_p05 <- list()

#create a data frame to store thresholds
threshold_summary <- data.frame(
  GWAS_POP = character(),
  DB_POP = character(),
  CELL_TYPE = character(),
  N_GENES = integer(),
  THRESHOLD = numeric(),
  stringsAsFactors = FALSE
)

#parse through files
for(file in files) {
  
  #read in file
  f <- fread(file)
  
  #unique identifier for this file (used as list key)
  model_name <- str_remove(basename(file), "\\.csv$")
  
  #assign variable (this will need to be changed as more ancestries/datasets are examined)
  gwas_pop <- "EUR"
  db_pop <- "EUR"
  cell_type <- args$cell_type
  n_genes <- nrow(f)
  bonferroni_threshold <- 0.05 / n_genes

  cat("Processing:", model_name, "\n")
  cat("  GWAS Pop:", gwas_pop, "| DB Pop:", db_pop, "| Cell Type:", cell_type)
  cat("  Genes:", n_genes, "| Threshold:", format(bonferroni_threshold, scientific = TRUE), "\n\n")
  
  #add to threshold summary
  threshold_summary <- rbind(threshold_summary, 
                             data.frame(
                               GWAS_POP = gwas_pop,
                               DB_POP = db_pop,
                               CELL_TYPE = cell_type,
                               N_GENES = n_genes,
                               THRESHOLD = bonferroni_threshold
                             ))
  
  #filter by Bonferroni threshold
  filtered <- f %>%
    select(gene, gene_name, zscore, effect_size, pvalue, var_g, pred_perf_r2, 
           pred_perf_pval, pred_perf_qval, n_snps_used, n_snps_in_cov, n_snps_in_model) %>%
    filter(pvalue < bonferroni_threshold) %>%
    mutate(
      gwas_pop = gwas_pop,
      db_pop = db_pop,
      cell_type = cell_type
    )
  
  #filter by p < 0.05
  filtered_p05 <- f %>%
    select(gene, gene_name, zscore, effect_size, pvalue, var_g, pred_perf_r2, 
           pred_perf_pval, pred_perf_qval, n_snps_used, n_snps_in_cov, n_snps_in_model) %>%
    filter(pvalue < 0.05) %>%
    mutate(
      gwas_pop = gwas_pop,
      db_pop = db_pop,
      cell_type = cell_type
    )
  
  #add to results list
  if(nrow(filtered) > 0) {
    all_results[[model_name]] <- filtered
  }
  if(nrow(filtered_p05) > 0) {
    all_results_p05[[model_name]] <- filtered_p05
  }
}

#combine all results into one table
merged_results <- rbindlist(all_results, fill = TRUE)
merged_results_p05 <- rbindlist(all_results_p05, fill = TRUE)

#reorder columns for better readability
if(nrow(merged_results) > 0) {
  merged_results <- merged_results %>%
    select(gene, gene_name, gwas_pop, db_pop, cell_type,
           zscore, effect_size, pvalue, var_g, pred_perf_r2, 
           pred_perf_pval, pred_perf_qval, n_snps_used, n_snps_in_cov, n_snps_in_model)
}

if(nrow(merged_results_p05) > 0) {
  merged_results_p05 <- merged_results_p05 %>%
    select(gene, gene_name, gwas_pop, db_pop, cell_type,
           zscore, effect_size, pvalue, var_g, pred_perf_r2, 
           pred_perf_pval, pred_perf_qval, n_snps_used, n_snps_in_cov, n_snps_in_model)
}

#write merged output
output_file <- paste0("/home/ckrueger2/sc_pwas_in_ukb/filtered_output/merged_significant_results_bcorr_", args$cell_type, ".tsv")
output_file_p05 <- paste0("/home/ckrueger2/sc_pwas_in_ukb/filtered_output/merged_significant_results_p05_", args$cell_type, ".tsv")
threshold_file <- paste0("/home/ckrueger2/sc_pwas_in_ukb/filtered_output/bonferroni_thresholds_", cell_type, ".txt")

write.table(merged_results, output_file, row.names=FALSE, quote=FALSE, sep="\t")
write.table(merged_results_p05, output_file_p05, row.names=FALSE, quote=FALSE, sep="\t")
write.table(threshold_summary, threshold_file, row.names=FALSE, quote=FALSE, sep="\t")

#print summary statistics
cat("\n========== SUMMARY ==========\n")
cat("Total files processed:", length(files), "\n")
cat("Bonferroni significant genes:", nrow(merged_results), "\n")
cat("P < 0.05 genes:", nrow(merged_results_p05), "\n")
