#!/usr/bin/env Rscript
library(tidyverse)

# command line args -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
# output
in_csv     <- args[1]
out_tsv     <- args[2]

# read and process genome accessions --------------------------------------

gather <- read_csv(in_csv) %>%
  mutate(accession = str_extract(name, "GCA_\\d+\\.\\d")) %>%
  select(accession)
write_tsv(gather, out_tsv, col_names = F)