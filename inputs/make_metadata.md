On 08/01/2023, we searched the SRA using txid6943[Organism:noexp] AND ("biomol rna"[Properties] AND ("platform pacbio smrt"[Properties] OR "platform illumina"[Properties])). 
We exported the results to a text file using the "Send to" drop down menu, where we selected "File" and "Summary".
We edited the column names, removing special characaters (comma), replacing spaces with underscores, and converting all uppercase letters to lower case.
We also added the columns bioproject, publication_url, publication_doi, publication_title, excluded, oligo_dt_annotated, blood_meal_hour, sex, contam, pool_size, tissue, source

```
mamba install pysradb=2.1.0

cut -f6 inputs/metadata.tsv | sort | uniq | head -n 13 | xargs echo
# outputs: 
# SRP032795 SRP051699 SRP052078 SRP052091 SRP052106 SRP052108 SRP052114 SRP052123 SRP052145 SRP052154 SRP091404 SRP373454 SRP446981
pysradb metadata --saveto tmp.tsv --detailed SRP032795 SRP051699 SRP052078 SRP052091 SRP052106 SRP052108 SRP052114 SRP052123 SRP052145 SRP052154 SRP091404 SRP373454 SRP446981
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

By hand, fill in missing sample_title values
