library(data.table)
library(dplyr)
library(stringr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--cell_type", help="single-cell model cell type")
args <- parser$parse_args()

base_dir_logs <- paste0("/home/ckrueger2/logs/", args$cell_type)
base_dir_outputs <- paste0("/home/ckrueger2/sc_pwas_in_ukb/output/", args$cell_type, "_cell")

logs <- data.table(file = list.files(base_dir_logs, pattern = paste0(".*_", args$cell_type, "\\.(out|log)$")))
outputs <- data.table(file = list.files(base_dir_outputs, pattern = paste0(".*_", args$cell_type, "_predixcan_output\\.csv$")))
#output_oids <- fread("/home/ckrueger2/sc_pwas_in_ukb/oids.txt", header = FALSE)
#colnames(output_oids) <- "OlinkID"

head(logs)
head(outputs)
#head(output_oids)

cat("Log file pattern:", logs$file[1], "\n")
cat("Output file pattern:", outputs$file[1], "\n")

log_oids <- logs %>%
  mutate(
    OlinkID = str_extract(file, "OID\\d{5}(\\.\\d+)?")) %>%
    filter(!is.na(OlinkID)) %>%
    select(OlinkID)

output_oids <- outputs %>%
  mutate(
    OlinkID = str_extract(file, "OID\\d{5}(\\.\\d+)?")) %>%
    filter(!is.na(OlinkID)) %>%
    select(OlinkID)

missing <- log_oids[!log_oids$OlinkID %in% output_oids$OlinkID]
#missing <- output_oids[!output_oids$OlinkID %in% log_oids$OlinkID]

cat("Missing count:", nrow(missing), "\n")

head(missing)

write.table(missing, paste0("/home/ckrueger2/sc_pwas_in_ukb/oids_", args$cell_type, "_oom.txt"), row.names=FALSE, col.names=FALSE, sep="\n", quote=FALSE)
