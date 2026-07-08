library(dplyr)
library(data.table)
library(biomaRt)

#protein mapping file from UKB-PPP
mapping_file <- fread("~/sc_pwas_in_ukb/olink_protein_map_3k_v1.tsv") #from UKB-PPP synapse
mapping_filtered <- mapping_file %>%
  dplyr::select(OlinkID, UniProt, Assay, ensembl_id, chr, gene_start, gene_end)
write.table(mapping_filtered, "~/sc_pwas_in_ukb/mapped_proteins.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

#find genes in .db files (within ROCKS cluster)
#   for db in /home/ckrueger2/sc_pwas_in_ukb/models/*/predict_db_Model_training_filtered.db; do celltype=$(basename $(dirname "$db")); sqlite3 "$db" "SELECT gene FROM extra;" > ~/sc_pwas_in_ukb/gene_lists/${celltype}_genes.txt; done
#   scp /home/ckrueger2/sc_pwas_in_ukb/gene_lists/* claudia@compbio1.cs.luc.edu:/home/claudia/sc_gene_lists/

#isolate unique ENSGs
files <- list.files("~/sc_gene_lists", pattern = "_genes.txt$", full.names = TRUE)
genes <- unique(unlist(lapply(files, readLines)))
writeLines(sort(genes), "~/sc_pwas_in_ukb/all_genes_unique.txt")

#build biomart mapping
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37)

#query based on ENSG IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "transcription_start_site",
                 "transcript_start", "transcript_end", "ensembl_gene_id_version",
                 "ensembl_transcript_id_version", "external_gene_name"),
  filters = c("ensembl_gene_id", "transcript_gencode_basic"),
  values = list(genes, TRUE),
  mart = mart
)

#keep only results on reference chromosomes
results <- results[grepl("^[0-9XY]+$", results$chromosome_name), ]

#ensure transcription start site equals transcript start, reordering where necessary
results$transcription_start_site <- ifelse(
  results$transcription_start_site != results$transcript_start,
  results$transcript_start,
  results$transcription_start_site
)
results <- results[order(results$chromosome_name, results$transcription_start_site), ]

#keep only the largest transcript per gene
results <- results %>%
  mutate(transcript_length = transcript_end - transcript_start) %>%
  group_by(ensembl_gene_id) %>%
  slice_max(transcript_length, n = 1, with_ties = FALSE) %>%
  ungroup()

results_filtered <- results %>%
  group_by(ensembl_gene_id) %>%
  slice(which.max(transcript_length)) %>%
  ungroup()

#move ensembl gene id column to first position
results_filtered <- results_filtered %>%
  dplyr::select(ensembl_gene_id, everything())

#sort based on chromosome and transcription start site
chromosome_order <- c(as.character(1:22), "X", "Y", "MT")
results_filtered$chromosome_name <- factor(results_filtered$chromosome_name, levels = chromosome_order)
results_filtered <- results_filtered[order(results_filtered$chromosome_name, results_filtered$transcription_start_site), ]

#save
write.table(results_filtered, "~/sc_pwas_in_ukb/mapped_genes.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
not_in_results <- genes[!(genes %in% results_filtered$ensembl_gene_id)]
print(not_in_results)

#parse through outputs
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(fs)

mapped_genes <- fread("/home/claudia/sc_pwas_in_ukb/mapped_genes.txt") %>%
  dplyr::select(ensembl_gene_id, chromosome_name, gene_tss = transcription_start_site, gene_stop = transcript_end)

mapped_proteins <- read_table("/home/claudia/sc_pwas_in_ukb/mapped_proteins.txt") %>%
  rename(olink_id = OlinkID, uniprot = UniProt,
         protein_tss = gene_start, protein_stop = gene_end) %>%
  dplyr::select(olink_id, uniprot, Assay, chr, protein_tss, protein_stop)

csv_files <- dir_ls("/home/claudia/sc_output/", recurse = TRUE, glob = "*_predixcan_output.csv") #keep all output in one folder (can contain subfolders), add subfolder to path if fdrs calculated by cell type

predixcan_all <- map_dfr(csv_files, function(f) {
  oid <- str_extract(basename(f), "OID\\d+")
  cell_type <- path_file(path_dir(f))
  read_csv(f, show_col_types = FALSE) %>%
    mutate(olink_id = oid, cell_type = cell_type,
           ensembl_gene_id = str_remove(gene, "\\..*"))
})

final_table <- predixcan_all %>%
  left_join(mapped_genes, by = "ensembl_gene_id") %>%
  left_join(mapped_proteins, by = "olink_id") %>%
  dplyr::select(olink_id, uniprot, ensembl_id = ensembl_gene_id, gene_name, cds_name = Assay, 
                gene_chr = chromosome_name, protein_chr = chr, gene_tss, protein_tss, gene_stop, protein_stop,
                zscore, effect_size, pvalue, cell_type)

#annotate with cis/trans-acting, cis-same/different
annotated_table <- final_table %>%
  mutate(
    acting = if_else(gene_chr == protein_chr & abs(gene_tss - protein_tss) <= 1e6, "cis", "trans"),
    cis_category = case_when(
      acting == "trans" ~ NA_character_,
      gene_name == cds_name ~ "same",
      TRUE ~ "different"
    )
  )

#calculate qvalues
library(qvalue)

all_pvals <- annotated_table$pvalue
qobj <- qvalue(all_pvals)
all_qvals <- qobj$qvalues
annotated_table$qvalue <- all_qvals
annotated_table <- annotated_table %>%
  relocate(qvalue, .after = pvalue)

write.table(annotated_table, "~/sc_pwas_in_ukb/sc_predixcan_results_annotated.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

#filter to FDR < 0.05
fdr_sig <- annotated_table %>%
  dplyr::filter(qvalue < 0.05) %>%
  dplyr::select(olink_id, uniprot, ensembl_id, gene_name, cds_name, zscore, effect_size, pvalue, qvalue, cell_type, acting, cis_category)

write.table(fdr_sig, "~/sc_pwas_in_ukb/sc_predixcan_results_fdr.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
