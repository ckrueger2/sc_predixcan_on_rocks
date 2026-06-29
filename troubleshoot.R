library(dplyr)
library(data.table)
library(stringr)

file_names <- list.files("/home/claudia/ukbio/euro_amer/")
file_oids <- str_extract(file_names, "OID\\d+")

write.table(file_oids, "/home/claudia/sc_pwas_in_ukb/oids.txt", quote=FALSE, sep="\n", col.names=FALSE, row.names=FALSE)
