On 08/01/2023, we searched the SRA using txid6943[Organism:noexp] AND ("biomol rna"[Properties] AND ("platform pacbio smrt"[Properties] OR "platform illumina"[Properties])). 
We exported the results to a text file using the "Send to" drop down menu, where we selected "File" and "Summary".
We edited the column names, removing special characaters (comma), replacing spaces with underscores, and converting all uppercase letters to lower case.
We also added the columns bioproject, publication_url, publication_doi, publication_title, excluded, oligo_dt_annotated, blood_meal_hour, sex, contam, pool_size, tissue, source

```
mamba install pysradb=2.1.0
# the cut command outputs the unique SRP accessions in inputs/metadata.tsv in the format that pysradb expects them in
pysradb metadata --saveto tmp.tsv --detailed $(cut -f6 inputs/metadata.tsv | sort | uniq | head -n 13 | xargs echo)
```

```
R
library(tidyverse)
m1 <- read_tsv("inputs/metadata.tsv")
m2 <- read_tsv("tmp.tsv")
metadata <- left_join(m1, m2, by = c('experiment_accession', 'experiment_title', 
                                     'organism_name', 'instrument', 'study_accession', 
                                     'study_title', 'sample_accession', 'sample_title',
                                     'total_spots', 'library_name', 'library_strategy', 
                                     'library_source', 'library_selection'))
write_tsv(metadata, "inputs/metadata.tsv")
```

By hand, fill in missing sample_title values.

Also created assembly_group by hand.
I separated samples by:
* different chemistries (Illumina, PacBio)
* read type (paired end,  single end)
* different tissues/cell lines (whole, salivary gland, gut, cell lines)
* sex (male vs. female) 
