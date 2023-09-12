#!/usr/bin/env Rscript
library(tidyverse)

# command line args -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

in_blast          <- args[1]
in_fai            <- args[2]
out_clean         <- args[3]
out_endosymbionts <- args[4]

# parse blast results -----------------------------------------------------

blast <- read_tsv(in_blast, col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                                          "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  # filter to the best BLAST hit
  group_by(qseqid) %>%
  slice_max(bitscore) %>%
  slice_head(n = 1)

fai <- read_tsv(in_fai, col_names = c("qseqid", "transcript_length", "offset", 
                                      "linebases", "linewidth", "qualoffset"))

blast <- left_join(blast, fai, by = c("qseqid")) %>%
  mutate(hit_fraction = length/transcript_length)

blast_filtered <- blast %>%
  filter(hit_fraction > 0.1) %>%
  filter(length > 100) %>%
  filter(pident >= 80)

# filter to contig names that are not in the contaminant list
keep <- fai %>%
  filter(!qseqid %in% blast_filtered$qseqid) %>%
  select(qseqid)

endosymbiont <- blast_filtered %>%
  filter(sseqid %in% c("CP021379.1")) # accession for coxiella-like endosymbiont contig

# write out the clean contig names
write_tsv(keep, out_clean, col_names = F)

# write out the endosymbiont contig names
write_tsv(endosymbiont, out_endosymbionts, col_names = F)