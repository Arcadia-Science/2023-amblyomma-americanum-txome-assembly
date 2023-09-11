import pandas as pd

# read in metadata file
metadata_all = pd.read_csv("../../inputs/metadata.tsv", sep = "\t").set_index("run_accession", drop = False)
# filter out samples that should be excluded (library prep was weird)
metadata_all = metadata_all[metadata_all['excluded'] == "keep"]
# select columns that we need metadata from for wildcards and other places in the workflow
metadata_filt = metadata_all[["library_name", "assembly_group", "library_layout", "instrument"]]
# separate the isoseq data bc it will be treated separately
metadata_illumina = metadata_filt[metadata_filt["instrument"] != "Sequel II"].drop_duplicates()
metadata_isoseq = metadata_all[metadata_filt["instrument"] == "Sequel II"]
# set the index to library name to allow dictionary-like lookups from the metadata tables with param lambda functions
metadata_illumina = metadata_illumina.set_index("library_name", drop = False)
# set the index to assembly group to allow dictionary-like lookups from the metadata tables with param lambda functions
metadata_illumina2 = metadata_illumina[["assembly_group", "library_layout"]].drop_duplicates()
metadata_illumina2 = metadata_illumina2.set_index("assembly_group", drop = False)

# use metadata tables to create global variables
# extract SRA accessions to a variable, which we'll use to download the raw data
RUN_ACCESSIONS = metadata_all["run_accession"].unique().tolist()
# extract library names, which we'll use to control the quality control portion of the workflow
# some libraries are split between multiple SRA accessions
ILLUMINA_LIB_NAMES = metadata_illumina["library_name"].unique().tolist()
# extract assembly groups, which we'll use to control the asesmbly portion of the workflow
ASSEMBLY_GROUPS = metadata_illumina["assembly_group"].unique().tolist()

# extract isoseq library names
ISOSEQ_LIB_NAMES = metadata_isoseq['library_name'].unique().tolist()
ISOSEQ_RUN_ACCESSIONS = metadata_isoseq['run_accession'].unique().tolist()

# set the short read assemblers
# ASSEMBLERS = ["trinity", "rnaspades"]
ASSEMBLERS = ["rnaspades", "trinity"]

READS = ['R1', 'R2']

rule all:
    input: expand("outputs/salmon/{assembly_group}_quant/quant.sf", assembly_group = ASSEMBLY_GROUPS)

rule rename_contigs:
    """
    we should talk to austin and figure out what his requirements are for transcript FASTA headers for noveltree
    """
    input: "outputs/assembly/{assembler}/{assembly_group}_{assembler}.fa"
    output: "outputs/assembly/renamed/{assembly_group}_{assembler}_renamed.fa"
    conda: "envs/bbmap.yml"
    shell:'''
    bbrename.sh in={input} out={output} prefix={wildcards.assembly_group} addprefix=t
    '''

rule merge_txomes_all:
    input: 
        expand("outputs/assembly/renamed/{assembly_group}_{assembler}_renamed.fa", assembly_group = ASSEMBLY_GROUPS, assembler = ASSEMBLERS),
        expand("outputs/assembly/isoseq/{isoseq_lib_name}.fa", isoseq_lib_name = ISOSEQ_LIB_NAMES)
    output: "outputs/assembly/merged/merged.fa"
    shell:'''
    cat {input} > {output}
    '''

rule deduplicate_merged_txomes_with_mmseqs:
    '''
    remove perfect duplicates.
    emulates cd-hit: https://github.com/soedinglab/MMseqs2/issues/601
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
    grab fasta header names and remove the prefix ">"
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

rule orthofinder:
    """
    orthofuser runs this step with orthofuser.py. 
    Orthofuser.py is a modified version of orthofinder that accepted nucleotides before orthofinder was set up to do that.
    since orthofinder will now run on nucleotdies on its own, we can use that directly
    * -d:  input is DNA sequence
    * -f:  input folder
    * -I:  MCL parameter
    * -og: stop after inferring orthogroups
    * -t:  number of parallel sequence search threads [default = 10]
    * -a:  number of parallel analysis threads
    """
    input: 
        illumina = expand("outputs/assembly/filtered_size/{assembly_group}_{assembler}_filtered.fa", assembly_group = ASSEMBLY_GROUPS, assembler = ASSEMBLERS),
        isoseq = expand("outputs/assembly/isoseq/{isoseq_lib_name}.fa", isoseq_lib_name = ISOSEQ_LIB_NAMES)
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

# TODO: add a rule to keep anything smaller than 75 and then come up with a strategy to work with this set (search for peptides, etc).

rule merge_by_assembly_group_for_transrate_round1:
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
        #reads=expand("outputs/fastp_separated_reads_combined/{read}.fq.gz", read = READS)
        reads=expand("outputs/assembly_group_separated_reads/{{assembly_group}}_{read}.fq.gz", read = READS)
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
        fa = "outputs/assembly/isoseq/{isoseq_lib_name}.fa",
        db = "inputs/databases/swissprot.dmnd"
    output: "outputs/orthofuser/diamond/{isoseq_lib_name}.diamond.txt"
    conda: "envs/diamond.yml"
    threads: 28
    shell:'''
    diamond blastx --quiet -p {threads} -e 1e-8 --top 0.1 -q {input.fa} -d {input.db} -o {output}
    '''

rule parse_diamond_gene_annotations_for_missed_transcripts: 
    input:
        orthomerged = "outputs/orthofuser/diamond/orthomerged.diamond.txt",
        raw = expand("outputs/orthofuser/diamond/{assembly_group}_{assembler}.diamond.txt", assembly_group = ASSEMBLY_GROUPS, assembler = ASSEMBLERS),
        isoseq = expand("outputs/orthofuser/diamond/{isoseq_lib_name}.diamond.txt", isoseq_lib_name = ISOSEQ_LIB_NAMES)
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

# consider another round of transrate using the diginorm'd full read set on the first "orthofuser_final.fa" transcriptome

rule salmon_index:
    input: "outputs/orthofuser/orthofuser_final.fa"
    output: "outputs/salmon/orthofuser_final_index/info.json"
    threads: 8
    params: indexdir = "outputs/salmon/orthofuser_final_index/"
    conda: "envs/salmon.yml"
    shell:'''
    salmon index -p {threads} -t {input} -i {params.indexdir} -k 31
    '''

rule salmon_quant:
    input:
        index = "outputs/salmon/orthofuser_final_index/info.json",
        reads=expand("outputs/assembly_group_separated_reads/{{assembly_group}}_{read}.fq.gz", read = READS)
    output: "outputs/salmon/{assembly_group}_quant/quant.sf"
    params: 
        indexdir = "outputs/salmon/orthofuser_final_index/",
        outdir = lambda wildcards: "outputs/salmon/" + wildcards.assembly_group + "_quant" 
    conda: "envs/salmon.yml"
    threads: 4
    shell:'''
    salmon quant -i {params.indexdir} -l A -1 {input.reads[0]} -2 {input.reads[1]} -o {params.outdir} --dumpEq --writeOrphanLinks -p {threads} 
    '''

