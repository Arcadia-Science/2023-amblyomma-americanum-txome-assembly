import pandas as pd

# read in metadata file
metadata_all = pd.read_csv("inputs/metadata.tsv", sep = "\t").set_index("run_accession", drop = False)
# filter out samples that should be excluded (library prep was weird)
metadata_all = metadata_all[metadata_all['excluded'] == "keep"]
# select columns that we need metadata from for wildcards and other places in the workflow
metadata_filt = metadata_all[["library_name", "assembly_group", "library_layout", "instrument"]]
# separate the isoseq data bc it will be treated separately
metadata_illumina = metadata_filt[metadata_filt["instrument"] != "Sequel II"].drop_duplicates()
metadata_isoseq = metadata_filt[metadata_filt["instrument"] == "Sequel II"]
# set the index to library name to allow dictionary-like lookups from the metadata tables with param lambda functions
metadata_illumina = metadata_illumina.set_index("library_name", drop = False)
# set the index to assembly group to allow dictionary-like lookups from the metadata tables with param lambda functions
metadata_illumina2 = metadata_illumina.set_index("assembly_group", drop = False)

# use metadata tables to create global variables
# extract SRA accessions to a variable, which we'll use to download the raw data
RUN_ACCESSIONS = metadata_all["run_accession"].unique().tolist()
# extract library names, which we'll use to control the quality control portion of the workflow
# some libraries are split between multiple SRA accessions
ILLUMINA_LIB_NAMES = metadata_illumina["library_name"].unique().tolist()
# extract assembly groups, which we'll use to control the asesmbly portion of the workflow
ASSEMBLY_GROUPS = metadata_illumina["assembly_group"].unique().tolist()

rule all:
    input: 
        expand("outputs/assembly/trinity/{assembly_group}_trinity.fa", assembly_group = ASSEMBLY_GROUPS),
        expand("outputs/assembly/rnaspades/{assembly_group}_rnaspades_hard_filtered_transcripts.fa", assembly_group = ASSEMBLY_GROUPS)

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
    output: expand("outputs/raw_combined/{illumina_lib_name}.fq.gz", illumina_lib_name = ILLUMINA_LIB_NAMES)
    params: 
        indir = "inputs/raw/",
        outdir = "outputs/raw_combined/"
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
    input: "outputs/raw_combined/{illumina_lib_name}.fq.gz"
    output:
        json = "outputs/fastp/{illumina_lib_name}.json",
        html = "outputs/fastp/{illumina_lib_name}.html",
        fq = "outputs/fastp/{illumina_lib_name}.fq.gz"
    conda: "envs/fastp.yml"
    threads: 2
    params: liblayout = lambda wildcards: metadata_illumina.loc[wildcards.illumina_lib_name, "library_layout"]
    shell:'''
    if [ "{params.liblayout}" == "PAIRED" ]; then
        fastp -i {input} --thread {threads} --trim_poly_x --qualified_quality_phred 2 --json {output.json} --html {output.html} --report_title {wildcards.illumina_lib_name} --interleaved_in --stdout | gzip > {output.fq}
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        fastp -i {input} --thread {threads} --trim_poly_x --qualified_quality_phred 2 --json {output.json} --html {output.html} --report_title {wildcards.illumina_lib_name} --stdout | gzip > {output.fq}
    fi
    '''

rule khmer_kmer_trim_and_normalization:
    """
    K-mer trim and diginorm (digital normalization, or just normalization) according to the eelpond protocol/elvers.
    The oyster river protocol also supports removal of erroneous k-mers through similar methods.
    """
    input: "outputs/fastp/{illumina_lib_name}.fq.gz"
    output: "outputs/khmer/{illumina_lib_name}.fq.gz"
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
    input: expand("outputs/khmer/{illumina_lib_name}.fq.gz", illumina_lib_name = ILLUMINA_LIB_NAMES)
    output: expand("outputs/assembly_group_interleaved_reads/{assembly_group}.fq.gz", assembly_group = ASSEMBLY_GROUPS)
    params:
        indir = "outputs/khmer/",
        outdir = "outputs/assembly_group_interleaved_reads"
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
    input: "outputs/assembly_group_interleaved_reads/{assembly_group}.fq.gz"
    output: 
        r1="outputs/assembly_group_separated_reads/{assembly_group}_R1.fq.gz",
        r2="outputs/assembly_group_separated_reads/{assembly_group}_R2.fq.gz"
    conda: "envs/bbmap.yml"
    params: liblayout = lambda wildcards: metadata_illumina2.loc[wildcards.assembly_group, "library_layout"]
    shell:'''
    if [ "{params.liblayout}" == "PAIRED" ]; then
        repair.sh in={input} out={output.r1} out2={output.r2} repair=t overwrite=true
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        cp {input} {output.r1}
        touch {output.r2}
    fi
    '''

rule trinity_assemble:
    input:
        r1="outputs/assembly_group_separated_reads/{assembly_group}_R1.fq.gz",
        r2="outputs/assembly_group_separated_reads/{assembly_group}_R2.fq.gz"
    output: "outputs/assembly/trinity/{assembly_group}_trinity.fa"
    conda: "envs/trinity.yml"
    threads: 7
    params: 
        liblayout = lambda wildcards: metadata_illumina2.loc[wildcards.assembly_group, "library_layout"],
        outdir = lambda wildcards: "outputs/assembly/trinity_tmp/" + wildcards.assembly_group + "_Trinity" 
    shell:'''
    if [ "{params.liblayout}" == "PAIRED" ]; then
        Trinity --left {input.r1} --right {input.r2} --seqType fq --CPU {threads} --max_memory 16G --output {params.outdir}
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        Trinity --single {input.r1} --seqType fq --CPU {threads} --max_memory 16G --output {params.outdir} 
    fi
    mv {params.outdir}.Trinity.fasta {output}
    '''

rule rnaspades_assemble:
    input:
        r1="outputs/assembly_group_separated_reads/{assembly_group}_R1.fq.gz",
        r2="outputs/assembly_group_separated_reads/{assembly_group}_R2.fq.gz"
    output: 
        hard = "outputs/assembly/rnaspades/{assembly_group}_rnaspades_hard_filtered_transcripts.fa",
        soft = "outputs/assembly/rnaspades/{assembly_group}_rnaspades_soft_filtered_transcripts.fa",
    conda: "envs/spades.yml"
    threads: 4
    params: 
        liblayout = lambda wildcards: metadata_illumina2.loc[wildcards.assembly_group, "library_layout"],
        outdir = lambda wildcards: "outputs/assembly/rnaspades_tmp/" + wildcards.assembly_group 
    shell:'''
    if [ "{params.liblayout}" == "PAIRED" ]; then
        rnaspades.py -1 {input.r1} -2 {input.r2} -o {params.outdir} -t {threads}
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        rnaspades.py -s {input.r1} -o {params.outdir} -t {threads}
    fi
    mv {params.outdir}/hard_filtered_transcripts.fasta {output.hard}
    mv {params.outdir}/soft_filtered_transcripts.fasta {output.soft}
    '''
