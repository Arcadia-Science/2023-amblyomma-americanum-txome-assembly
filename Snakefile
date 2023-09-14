import pandas as pd

# read in metadata file
metadata_all = pd.read_csv("inputs/metadata.tsv", sep = "\t").set_index("run_accession", drop = False)
# filter out samples that should be excluded (library prep was weird)
metadata_all = metadata_all[metadata_all['excluded'] == "keep"]
# select columns that we need metadata from for wildcards and other places in the workflow
metadata_filt = metadata_all[["library_name", "assembly_group", "library_layout", "instrument"]]
# separate the isoseq data bc it will be treated separately
metadata_illumina = metadata_filt[metadata_filt["instrument"] != "Sequel II"].drop_duplicates()
metadata_isoseq = metadata_all[metadata_filt["instrument"] == "Sequel II"]
# set the index to library name to allow dictionary-like lookups from the metadata tables
metadata_illumina = metadata_illumina.set_index("library_name", drop = False)

# use metadata tables to create global variables
# extract SRA accessions to a variable, which we'll use to download the raw data
RUN_ACCESSIONS = metadata_all["run_accession"].unique().tolist()
# extract library names, which we'll use to control the quality control portion of the workflow
# some libraries are split between multiple SRA accessions
ILLUMINA_LIB_NAMES = metadata_illumina["library_name"].unique().tolist()
# extract assembly groups, which we'll use to control the assembly portion of the workflow
ASSEMBLY_GROUPS = metadata_illumina["assembly_group"].unique().tolist()

# extract isoseq library names
ISOSEQ_LIB_NAMES = metadata_isoseq['library_name'].unique().tolist()
ISOSEQ_RUN_ACCESSIONS = metadata_isoseq['run_accession'].unique().tolist()

# set the short read assemblers
ASSEMBLERS = ["rnaspades", "trinity"]

READS = ['R1', 'R2']

# set contam screen params
LINEAGES = ['bacteria', 'archaea', 'protozoa', 'fungi', 'vertebrate_mammalian', 'vertebrate_other', 'plant']
KSIZES = [51]

rule all:
    input: 
        expand("outputs/evaluation/salmon/{assembly_group}_quant/quant.sf", assembly_group = ASSEMBLY_GROUPS), 
        "outputs/decontamination/orthofuser_final_endosymbiont.fa",
        "outputs/annotation/transdecoder/orthofuser_final_clean.fa.transdecoder.cds"

######################################
# Download short & long read data
######################################

rule download_fastq_files:
    output: temp("inputs/raw/{run_accession}.fq.gz")
    conda: "envs/sratools.yml"
    params: outdir = "inputs/raw/"
    shell:'''
    fasterq-dump --split-spot -Z {wildcards.run_accession} | gzip > {output}
    ''' 

######################################
## Process & assemble illumina files
######################################

rule combine_by_library_name:
    """
    Some of the input sequences have multiple run accessions (SRR*) associated with a single library.
    This rule combines those run accessions into one file by ill_lib_name (illumina library name).
    Since it uses the metadata_illumina file to do this, as well as expanding over the runacc wild card in the input,
    the output wildcard illumina_lib_name will only contain illumina library names.
    """
    input: expand("inputs/raw/{run_accession}.fq.gz", run_accession = RUN_ACCESSIONS)
    output: temp(expand("outputs/read_qc/raw_combined/{illumina_lib_name}.fq.gz", illumina_lib_name = ILLUMINA_LIB_NAMES))
    params: 
        indir = "inputs/raw/",
        outdir = "outputs/read_qc/raw_combined/"
    run:
        # create a dictionary that has library names as keys and run accessions as values
        tmp_illumina = metadata_all[metadata_all["instrument"] != "Sequel II"].drop_duplicates()
        tmp = tmp_illumina[["library_name", "run_accession"]]
        libdict = {}
        for group, d in tmp.groupby('library_name'):
            libdict[group] = d['run_accession'].values.tolist()

        # format the values to be a string of file paths, each separated by a space
        for library_name, run_accessions in libdict.items():
            library_paths = []
            for run_accession in run_accessions:
                path = params.indir + run_accession + ".fq.gz"
                library_paths.append(path)
                shell_drop_in = " ".join(library_paths)
            shell("cat {shell_drop_in} > {params.outdir}/{library_name}.fq.gz")

rule fastp:
    """
    We set the quality trimming parameters used by the Oyster River Protocol for de novo transcriptomics:
    - quality trim with a phred score of 2
    - trim the polyA tails
    - adapter trim
    """
    input: "outputs/read_qc/raw_combined/{illumina_lib_name}.fq.gz"
    output:
        json = "outputs/read_qc/fastp/{illumina_lib_name}.json",
        html = "outputs/read_qc/fastp/{illumina_lib_name}.html",
        fq = "outputs/read_qc/fastp/{illumina_lib_name}.fq.gz"
    conda: "envs/fastp.yml"
    threads: 2
    shell:'''
    fastp -i {input} --thread {threads} --trim_poly_x --qualified_quality_phred 2 --json {output.json} --html {output.html} --report_title {wildcards.illumina_lib_name} --interleaved_in --stdout | gzip > {output.fq}
    '''

rule khmer_kmer_trim_and_normalization:
    """
    K-mer trim and diginorm (digital normalization, or just normalization) according to the eelpond protocol/elvers.
    The oyster river protocol also supports removal of erroneous k-mers through similar methods.
    """
    input: "outputs/read_qc/fastp/{illumina_lib_name}.fq.gz"
    output: "outputs/read_qc/khmer/{illumina_lib_name}.fq.gz"
    conda: "envs/khmer.yml"
    shell:'''
    trim-low-abund.py -V -k 20 -Z 18 -C 2 -o {output} -M 4e9 --diginorm --diginorm-coverage=20 --gzip {input}
    '''

rule combine_by_assembly_group:
    """
    Assembly is a balancing act for de novo transcriptomics.
    Read depth must be sufficient to maximize coverage of rarely expressed transcripts,
    while isoform and SNP variation should be decreased as much as possible.
    This rule combines RNA-seq samples by pre-determined assembly groups (see metadata['assembly_group']) selected to balance the two above constraints.
    """
    input: expand("outputs/read_qc/khmer/{illumina_lib_name}.fq.gz", illumina_lib_name = ILLUMINA_LIB_NAMES)
    output: expand("outputs/read_qc/assembly_group_interleaved_reads/{assembly_group}.fq.gz", assembly_group = ASSEMBLY_GROUPS)
    params:
        indir = "outputs/read_qc/khmer/",
        outdir = "outputs/read_qc/assembly_group_interleaved_reads"
    run:
        # create a dictionary that has assembly groups as keys and library names as values
        tmp = metadata_illumina[["assembly_group", "library_name"]]
        assembly_group_dict = {}
        for group, d in tmp.groupby('assembly_group'):
            assembly_group_dict[group] = d['library_name'].values.tolist()

        # format the values to be a string of file paths, each separated by a space
        for assembly_group, library_names in assembly_group_dict.items():
            assembly_group_paths = []
            for library_name in library_names:
                path = params.indir + library_name + ".fq.gz"
                assembly_group_paths.append(path)
            
            shell_drop_in = " ".join(assembly_group_paths)
            shell("cat {shell_drop_in} > {params.outdir}/{assembly_group}.fq.gz")


rule split_paired_end_reads:
    """
    The trinity transcriptome assembler don't take interleaved reads as input.
    This rule separates reads into forward (R1) and reverse (R2) pairs.
    For single end reads, it touches the R2 file, as there is no information to separate.
    """
    input: "outputs/read_qc/assembly_group_interleaved_reads/{assembly_group}.fq.gz"
    output: 
        r1="outputs/read_qc/assembly_group_separated_reads/{assembly_group}_R1.fq.gz",
        r2="outputs/read_qc/assembly_group_separated_reads/{assembly_group}_R2.fq.gz"
    conda: "envs/bbmap.yml"
    shell:'''
    repair.sh in={input} out={output.r1} out2={output.r2} repair=t overwrite=true
    '''

rule trinity_assemble:
    input:
        r1="outputs/read_qc/assembly_group_separated_reads/{assembly_group}_R1.fq.gz",
        r2="outputs/read_qc/assembly_group_separated_reads/{assembly_group}_R2.fq.gz"
    output: "outputs/assembly/trinity/{assembly_group}_trinity.fa"
    conda: "envs/trinity.yml"
    threads: 28
    params: outdir = lambda wildcards: "outputs/assembly/trinity_tmp/" + wildcards.assembly_group + "_Trinity" 
    shell:'''
    Trinity --left {input.r1} --right {input.r2} --seqType fq --CPU {threads} --max_memory 100G --output {params.outdir} --full_cleanup
    mv {params.outdir}.Trinity.fasta {output}
    '''

rule rnaspades_assemble:
    input:
        r1="outputs/read_qc/assembly_group_separated_reads/{assembly_group}_R1.fq.gz",
        r2="outputs/read_qc/assembly_group_separated_reads/{assembly_group}_R2.fq.gz"
    output: 
        hard = "outputs/assembly/rnaspades/{assembly_group}_rnaspades_hard_filtered_transcripts.fa",
        soft = "outputs/assembly/rnaspades/{assembly_group}_rnaspades.fa"
    conda: "envs/spades.yml"
    threads: 4
    params: outdir = lambda wildcards: "outputs/assembly/rnaspades_tmp/" + wildcards.assembly_group 
    shell:'''
    rnaspades.py -1 {input.r1} -2 {input.r2} -o {params.outdir} -t {threads}
    mv {params.outdir}/hard_filtered_transcripts.fasta {output.hard}
    mv {params.outdir}/soft_filtered_transcripts.fasta {output.soft}
    '''

rule split_paired_end_reads_fastp:
    input: fq = "outputs/read_qc/fastp/{illumina_lib_name}.fq.gz"
    output:
        r1="outputs/read_qc/fastp_separated_reads/{illumina_lib_name}_R1.fq.gz",
        r2="outputs/read_qc/fastp_separated_reads/{illumina_lib_name}_R2.fq.gz"
    conda: "envs/bbmap.yml"
    shell:'''
    repair.sh in={input} out={output.r1} out2={output.r2} repair=t overwrite=true
    '''

######################################
## Process & assemble isoseq files
######################################

# The A americanum isoseq data on the SRA has already been processed.
# In this case, the sample was processed by UC Berkeley's computational core using the PacBio endorsed workflow.
#
# This probably looks something like this:
# 1. Generate CCS consensuses from raw isoseq subreads (bam file) (PBCCS)
# 2. Remove primer sequences from consensuses (LIMA)
# 3. Detect and remove chimeric reads (ISOSEQ3 REFINE)
# 4. Convert bam file into fasta file (BAMTOOLS CONVERT)
# 5. Select reads with a polyA tail and trim it (GSTAMA_POLYACLEANUP)
#
# Since these steps have already been completed, the FASTQ file we're working with here already represents a non-redundant set of the longest transcripts that could be derived from the raw data.
# We therefore only need to transform it into a FASTA file in order to include it in this analysis.
#
# In the future, if we need to do isoseq data processing, the first half of the nf-core/isoseq workflow has this above pipeline implemented.

rule convert_isoseq_fastq_to_fasta:
    input: expand("inputs/raw/{isoseq_run_accession}.fq.gz", isoseq_run_accession = ISOSEQ_RUN_ACCESSIONS)
    output: "outputs/assembly/isoseq/{isoseq_lib_name}_isoseq.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk seq -a {input} > {output}
    '''

##################################################
## Merge and deduplicate transcriptome assemblies
##################################################

rule rename_contigs:
    """
    prepends assembly group to each contig fasta header. 
    makes contig names that are duplicated between assemblies unique again.
    we should talk to austin and figure out what his requirements are for transcript FASTA headers for noveltree
    """
    input: "outputs/assembly/{assembler}/{assembly_group}_{assembler}.fa"
    output: "outputs/assembly/renamed/{assembly_group}_{assembler}_renamed.fa"
    conda: "envs/bbmap.yml"
    shell:'''
    bbrename.sh in={input} out={output} prefix={wildcards.assembly_group} addprefix=t
    '''

rule merge_txomes_all:
    '''
    combine all assembled contigs into one file
    '''
    input:
        expand("outputs/assembly/renamed/{assembly_group}_{assembler}_renamed.fa", assembly_group = ASSEMBLY_GROUPS, assembler = ASSEMBLERS),
        expand("outputs/assembly/isoseq/{isoseq_lib_name}_isoseq.fa", isoseq_lib_name = ISOSEQ_LIB_NAMES)
    output: "outputs/assembly/merged/merged.fa"
    shell:'''
    cat {input} > {output}
    '''

rule deduplicate_merged_txomes_with_mmseqs:
    '''
    Remove perfect duplicates (emulates cd-hit: https://github.com/soedinglab/MMseqs2/issues/601)
    The goal of this step is to remove fragments in one transcriptome that are present in a different transcriptome as longer contigs.
    '''
    input: "outputs/assembly/merged/merged.fa"
    output: "outputs/assembly/merged/merged_rep_seq.fasta"
    params: outprefix="outputs/assembly/merged/merged"
    conda: "envs/mmseqs2.yml"
    threads: 8
    shell:'''
    mkdir -p tmp
    mmseqs easy-cluster {input} {params.outprefix} tmp -c 0.97 --cov-mode 1 --min-seq-id 1.0 --exact-kmer-matching 1 --threads {threads}
    '''

rule grab_deduplicated_contig_names:
    """
    Grab fasta header names and remove the prefix ">" for deduplicated contigs
    """
    input: "outputs/assembly/merged/merged_rep_seq.fasta"
    output: "outputs/assembly/merged/merged_rep_seq_names.txt"
    shell:'''
    grep -e ">" {input} | awk 'sub(/^>/, "")' > {output}
    '''

rule filter_original_assemblies_by_deduplicated_names:
    """
    remove perfect duplicates from each original txome before running transrate or orthofinder.
    this step shares info between transcriptomes and is intended to remove small fragments that are part of larger transcripts in different assemblies and remove perfect duplicates between different assemblers.
    filtering is fast, while both transrate and orthofinder are slow.
    """
    input:
        fa="outputs/assembly/renamed/{assembly_group}_{assembler}_renamed.fa",
        lst="outputs/assembly/merged/merged_rep_seq_names.txt"
    output: "outputs/assembly/filtered_duplicates/{assembly_group}_{assembler}_filtered.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk subseq {input.fa} {input.lst} > {output}
    '''

rule filter_by_length:
    """
    Orthofuser filters to nucleotides greater than 200bp
    We'll use 75, since we're using 25 amino acids as our cut off (which syncs with what noveltree uses)
    """
    input: "outputs/assembly/filtered_duplicates/{assembly_group}_{assembler}_filtered.fa"
    output: "outputs/assembly/filtered_size/{assembly_group}_{assembler}_filtered.fa"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit seq -m 75 -o {output} {input}
    '''

# TODO: add a rule to keep anything smaller than 75 and then come up with a strategy to work with this set (search for peptides, etc).

rule orthofinder:
    """
    Using each assembly as an input "species," we run orthofinder to detect orthologous groups of transcripts.
    This is run on DNA sequences with an inflated MCL parameter to co-group more distant sequences.
    Orthofinder results (e.g. what is grouped together) might change if the input files change, so we run orthofinder one time to get clusters with consistent results and names.
    orthofinder parameters
    * -d:  input is DNA sequence
    * -f:  input folder
    * -I:  MCL parameter
    * -og: stop after inferring orthogroups
    * -t:  number of parallel sequence search threads [default = 10]
    * -a:  number of parallel analysis threads
    """
    input:
        illumina = expand("outputs/assembly/filtered_size/{assembly_group}_{assembler}_filtered.fa", assembly_group = ASSEMBLY_GROUPS, assembler = ASSEMBLERS),
        isoseq = expand("outputs/assembly/isoseq/{isoseq_lib_name}_isoseq.fa", isoseq_lib_name = ISOSEQ_LIB_NAMES)
    output: "outputs/orthofuser/orthofinder/Orthogroups/Orthogroups.txt"
    params:
        indir="outputs/assembly/filtered_size/",
        outdirtmp = "outputs/orthofuser/orthofinder_tmp"
    threads: 28
    conda: "envs/orthofinder.yml"
    shell:'''
    # put the isoseq data in the same folder as the illumina assemblies
    # in this case, since we only have one isoseq file, we can do a simple cp instead of a for loop
    cp {input.isoseq} {params.indir}
    orthofinder -d -I 12 -f {params.indir} -o {params.outdirtmp} -og -t {threads}
    cp {params.outdirtmp}/Results*/Orthogroups/Orthogroups.txt {output}
    '''

rule merge_by_assembly_group_for_transrate:
    input:
        expand("outputs/assembly/filtered_duplicates/{{assembly_group}}_{assembler}_filtered.fa", assembler = ASSEMBLERS),
    output: "outputs/assembly/merged/{assembly_group}_merged_filtered.fa"
    shell:'''
    cat {input} > {output}
    '''

rule transrate:
    """
    TransRate is a tool for reference-free quality assessment of de novo transcriptome assemblies.
    It uses evidence like read mapping rate and length to score each contig.
    orthofuser uses a modified version of transrate,
    but it only uses an updated version of salmon and other bugs, which shouldn't change the functionality
    https://github.com/macmanes-lab/Oyster_River_Protocol/issues/46

    Ideally, we would run transrate on the merged transcriptome and map with merged reads.
    For this transcriptome, the assembly is too complex (est. 1M contig limit) and there are too many reads, both of which cause transrate to fail.
    For this reason, we run transrate separately on each assembly group and then combine the scores to get the contig scores.
    Justification for using diginorm'd reads for mapping: https://github.com/blahah/transrate/issues/225
    """
    input:
        assembly="outputs/assembly/merged/{assembly_group}_merged_filtered.fa",
        reads=expand("outputs/read_qc/assembly_group_separated_reads/{{assembly_group}}_{read}.fq.gz", read = READS)
    output: "outputs/orthofuser/transrate_full/{assembly_group}_merged_filtered/contigs.csv"
    singularity: "docker://macmaneslab/orp:2.3.3"
    params: outdir= "outputs/orthofuser/transrate_full"
    threads: 28
    shell:'''
    transrate -o {params.outdir} -t {threads} -a {input.assembly} --left {input.reads[0]} --right {input.reads[1]}
    mv {params.outdir}/assemblies.csv {params.outdir}/{wildcards.assembly_group}_merged_filtered_assemblies.csv
    # cleanup big output files
    rm {params.outdir}/{wildcards.assembly_group}_merged_filtered/*bam
    '''

rule get_contig_name_w_highest_transrate_score_for_each_orthogroup:
    """
    This rule replaces a lot of the bash/awk/grep/thousands of file writing steps in orthofuser/ORP with a python script.
    It writes the name of the highest scoring transcript (contig name) from transrate for each orthogroup.
    It takes the highest scoring contig for each orthogroup from the multiple transrate runs.
    """
    input:
        orthogroups = "outputs/orthofuser/orthofinder/Orthogroups/Orthogroups.txt",
        transrate = expand("outputs/orthofuser/transrate_full/{assembly_group}_merged_filtered/contigs.csv", assembly_group = ASSEMBLY_GROUPS)
    output: "outputs/orthofuser/orthomerged/good.list"
    shell:'''
    python scripts/get_contig_name_w_highest_transrate_score_for_each_orthogroup.py {output} {input.orthogroups} {input.transrate}
    '''

rule filter_by_name:
    """
    keep only contigs that ended up in good.list
    """
    input:
        fa="outputs/assembly/merged/merged_rep_seq.fasta",
        lst="outputs/orthofuser/orthomerged/good.list"
    output: "outputs/orthofuser/orthomerged/orthomerged.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk subseq {input.fa} {input.lst} > {output}
    '''

rule download_diamond_database:
    output: "inputs/databases/uniprot_sprot.fasta.gz"
    shell:'''
    # downloaded UniProt Release 2023_03 on 20230818
    curl -JLo {output} http://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    '''

rule diamond_makedb:
    input: "inputs/databases/uniprot_sprot.fasta.gz"
    output: "inputs/databases/swissprot.dmnd"
    params: dbprefix = "inputs/databases/swissprot"
    conda: "envs/diamond.yml"
    shell:'''
    diamond makedb --in {input} --db {params.dbprefix}
    '''

rule run_diamond_on_orthomerged_txome:
    """
    get minimal annotations for the orthomerged transcriptome to see what genes are present.
    """
    input:
        fa ="outputs/orthofuser/orthomerged/orthomerged.fa",
        db = "inputs/databases/swissprot.dmnd"
    output: "outputs/orthofuser/diamond/orthomerged.diamond.txt"
    conda: "envs/diamond.yml"
    threads: 28
    shell:'''
    diamond blastx --quiet -p {threads} -e 1e-8 --top 0.1 -q {input.fa} -d {input.db} -o {output}
    '''

rule run_diamond_to_rescue_real_genes:
    """
    run diamond on the raw, unprocessed transcriptomes
    """
    input:
        fa = "outputs/assembly/filtered_duplicates/{assembly_group}_{assembler}_filtered.fa",
        db = "inputs/databases/swissprot.dmnd"
    output: "outputs/orthofuser/diamond/{assembly_group}_{assembler}.diamond.txt"
    conda: "envs/diamond.yml"
    threads: 28
    shell:'''
    diamond blastx --quiet -p {threads} -e 1e-8 --top 0.1 -q {input.fa} -d {input.db} -o {output}
    '''

rule run_diamond_to_rescue_real_genes_isoseq:
    """
    run diamond on the raw, unprocessed transcriptomes
    """
    input:
        fa = "outputs/assembly/isoseq/{isoseq_lib_name}_isoseq.fa",
        db = "inputs/databases/swissprot.dmnd"
    output: "outputs/orthofuser/diamond_isoseq/{isoseq_lib_name}_isoseq.diamond.txt"
    conda: "envs/diamond.yml"
    threads: 28
    shell:'''
    diamond blastx --quiet -p {threads} -e 1e-8 --top 0.1 -q {input.fa} -d {input.db} -o {output}
    '''

rule parse_diamond_gene_annotations_for_missed_transcripts:
    input:
        orthomerged = "outputs/orthofuser/diamond/orthomerged.diamond.txt",
        raw = expand("outputs/orthofuser/diamond/{assembly_group}_{assembler}.diamond.txt", assembly_group = ASSEMBLY_GROUPS, assembler = ASSEMBLERS),
        isoseq = expand("outputs/orthofuser/diamond_isoseq/{isoseq_lib_name}_isoseq.diamond.txt", isoseq_lib_name = ISOSEQ_LIB_NAMES)
    output: "outputs/orthofuser/newbies/newbies.txt"
    shell:'''
    python scripts/parse_diamond_gene_annotations_for_missed_transcripts.py {output} {input.orthomerged} {input.raw} {input.isoseq}
    '''

rule grab_missed_transcripts:
    input:
        fa="outputs/assembly/merged/merged_rep_seq.fasta",
        lst = "outputs/orthofuser/newbies/newbies.txt"
    output: "outputs/orthofuser/newbies/newbies.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk subseq {input.fa} {input.lst} > {output}
    '''

rule combine_orthomerged_with_missed_transcripts:
    input:
        orthomerged = "outputs/orthofuser/orthomerged/orthomerged.fa",
        newbies = "outputs/orthofuser/newbies/newbies.fa"
    output: "outputs/orthofuser/newbies/orthomerged.fa"
    shell:'''
    cat {input.orthomerged} {input.newbies} > {output}
    '''

rule cdhitest:
    """
    collapse at 0.98 identity, which should be removing nearly identical transcripts
    """
    input: "outputs/orthofuser/newbies/orthomerged.fa"
    output: "outputs/orthofuser/orthofuser_final.fa"
    conda: "envs/cd-hit.yml"
    threads: 30
    shell:'''
    cd-hit-est -M 5000 -T {threads} -c .98 -i {input} -o {output}
    '''

#######################################################################
## Evaluate the quality of the merged transcriptome
#######################################################################

rule salmon_index:
    input: "outputs/orthofuser/orthofuser_final.fa"
    output: "outputs/evaluation/salmon/orthofuser_final_index/info.json"
    threads: 8
    params: indexdir = "outputs/evaluation/salmon/orthofuser_final_index/"
    conda: "envs/salmon.yml"
    shell:'''
    salmon index -p {threads} -t {input} -i {params.indexdir} -k 31
    '''

rule salmon_quant:
    input:
        index = "outputs/evaluation/salmon/orthofuser_final_index/info.json",
        reads=expand("outputs/read_qc/assembly_group_separated_reads/{{assembly_group}}_{read}.fq.gz", read = READS)
    output: "outputs/evaluation/salmon/{assembly_group}_quant/quant.sf"
    params:
        indexdir = "outputs/evaluation/salmon/orthofuser_final_index/",
        outdir = lambda wildcards: "outputs/evaluation/salmon/" + wildcards.assembly_group + "_quant"
    conda: "envs/salmon.yml"
    threads: 4
    shell:'''
    salmon quant -i {params.indexdir} -l A -1 {input.reads[0]} -2 {input.reads[1]} -o {params.outdir} --dumpEq --writeOrphanLinks -p {threads}
    '''

## Screen for contamination ------------------------------------------

rule download_sourmash_databases_genbank:
    input: "inputs/sourmash_databases/sourmash-database-info.csv"
    output: "inputs/sourmash_databases/genbank-{lineage}-k{ksize}-scaled10k-cover.zip"
    run:
        sourmash_database_info = pd.read_csv(input[0])
        ksize = int(wildcards.ksize)
        lineage_df = sourmash_database_info.loc[(sourmash_database_info['lineage'] == wildcards.lineage) & (sourmash_database_info['ksize'] == ksize)]
        if lineage_df is None:
            raise TypeError("'None' value provided for lineage_df. Are you sure the sourmash database info csv was not empty?")

        osf_hash = lineage_df['osf_hash'].values[0] 
        shell("curl -JLo {output} https://osf.io/{osf_hash}/download")


rule sourmash_sketch:
    input: "outputs/orthofuser/orthofuser_final.fa"
    output: "outputs/decontamination/sourmash/orthofuser_final.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,abund,scaled=1000 -o {output} --name orthofuser_final {input}
    '''

rule sourmash_gather:
    input:
        sig="outputs/decontamination/sourmash/orthofuser_final.sig",
        db=expand("inputs/sourmash_databases/genbank-{lineage}-k{{ksize}}-scaled10k-cover.zip", lineage = LINEAGES),
    output: "outputs/decontamination/sourmash/orthofuser_final_k{ksize}_gather.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -k {wildcards.ksize} -o {output} {input.sig} {input.db}
    '''

rule parse_genome_accessions_from_sourmash_gather:
    input: expand("outputs/decontamination/sourmash/orthofuser_final_k{ksize}_gather.csv", ksize = KSIZES)
    output: "outputs/decontamination/orthofuser_final_contaminant_genome_accesions.tsv"
    conda: "envs/tidyverse.yml"
    shell:'''
    scripts/parse_genome_accessions_from_sourmash_gather_results.R {input} {output}
    '''

rule download_genome_accessions:
    input: "outputs/decontamination/orthofuser_final_contaminant_genome_accesions.tsv"
    output: "outputs/decontamination/contaminant_genomes/contaminant_genomes.fna"
    params: outdir = "outputs/decontamination/contaminant_genomes/"
    conda: "envs/ncbi-genome-download.yml"
    shell:'''
    ncbi-genome-download -s genbank --flat-output --formats fasta --output-folder {params.outdir} --assembly-accessions {input} all && cat {params.outdir}/*fna.gz | gunzip > {output}
    '''

rule make_contam_genome_blast_db:
    input: "outputs/decontamination/contaminant_genomes/contaminant_genomes.fna"
    output: "outputs/decontamination/contaminant_genomes_blastdb/contaminant_genomes_blastdb.not"
    params: dbprefix="outputs/decontamination/contaminant_genomes_blastdb/contaminant_genomes_blastdb"
    conda: "envs/blast.yml"
    shell:'''
    makeblastdb -in {input} -dbtype nucl -out {params.dbprefix}
    '''

rule run_blast_for_contam_screen:
    input: 
        fa="outputs/orthofuser/orthofuser_final.fa",
        db="outputs/decontamination/contaminant_genomes_blastdb/contaminant_genomes_blastdb.not"
    output: "outputs/decontamination/contaminant_blast/contaminant_genomes_blast.tsv"
    threads: 14
    params: dbprefix="outputs/decontamination/contaminant_genomes_blastdb/contaminant_genomes_blastdb"
    conda: "envs/blast.yml"
    shell:'''
    blastn -query {input.fa} -db {params.dbprefix} -outfmt 6 -out {output} -num_threads {threads} -max_target_seqs 5
    '''

rule samtools_faidx_transcriptome:
    input: "outputs/orthofuser/orthofuser_final.fa"
    output: "outputs/orthofuser/orthofuser_final.fa.fai"
    conda: "envs/samtools.yml"
    shell:'''
    samtools faidx {input}
    '''

rule parse_blast_for_contaminant_contigs:
    input:
        blast="outputs/decontamination/contaminant_blast/contaminant_genomes_blast.tsv",
        fai="outputs/orthofuser/orthofuser_final.fa.fai"
    output: 
        clean="outputs/decontamination/contaminant_blast/clean_contig_names.txt",
        endosymbiont="outputs/decontamination/contaminant_blast/endosymbiont_contig_names.txt"
    conda: "envs/tidyverse.yml"
    shell:'''
    scripts/parse_blast_for_contaminant_contigs.R {input.blast} {input.fai} {output.clean} {output.endosymbiont}
    '''

rule extract_clean_contigs:
    input: 
        fa="outputs/orthofuser/orthofuser_final.fa",
        clean="outputs/decontamination/contaminant_blast/clean_contig_names.txt",
    output: "outputs/decontamination/orthofuser_final_clean.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk subseq {input.fa} {input.clean} > {output}
    '''

rule extract_endosymbiont_contigs:
    input: 
        fa="outputs/orthofuser/orthofuser_final.fa",
        endosymbiont="outputs/decontamination/contaminant_blast/endosymbiont_contig_names.txt",
    output: "outputs/decontamination/orthofuser_final_endosymbiont.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk subseq {input.fa} {input.endosymbiont} > {output}
    '''

################################################
## ORF prediction
################################################

rule transdecoder_longorfs:
    input: "outputs/decontamination/orthofuser_final_clean.fa"
    output: "outputs/annotation/transdecoder/orthofuser_final_clean.fa.transdecoder_dir/longest_orfs.cds"
    params: outdir="outputs/annotation/transdecoder/"
    conda: "envs/transdecoder.yml"
    shell:'''
    TransDecoder.LongOrfs -t {input} --output_dir {params.outdir}
    '''

rule transdecoder_predict:
    input: 
        td="outputs/annotation/transdecoder/orthofuser_final_clean.fa.transdecoder_dir/longest_orfs.cds",
        fa="outputs/decontamination/orthofuser_final_clean.fa"
    output: "outputs/annotation/transdecoder/orthofuser_final_clean.fa.transdecoder.cds"
    conda: "envs/transdecoder.yml"
    params: outdir="outputs/annotation/transdecoder/"
    shell:'''
    TransDecoder.Predict -t {input.fa} --output_dir {params.outdir}
    '''
