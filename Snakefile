import pandas as pd

# read in metadata file
metadata_all = pd.read_csv("inputs/metadata.tsv", sep = "\t").set_index("run_accession", drop = False)
# select columns that we need metadata from for wildcards and other places in the workflow
metadata_filt = metadata_all[["library_name", "assembly_group", "library_layout", "instrument"]]
# separate the isoseq data bc it will be treated separately
metadata_illumina = metadata_filt[metadata_filt["instrument"] != "Sequel II"].drop_duplicates()
metadata_isoseq = metadata_filt[metadata_filt["instrument"] == "Sequel II"]
# set the index to library name to allow dictionary-like lookups from the metadata tables with param lambda functions
metadata_illumina = metadata_illumina.set_index("library_name", drop = False)

# use metadata tables to create global variables
# extract SRA accessions to a variable, which we'll use to download the raw data
RUN_ACCESSIONS = metadata_all["run_accession"].unique().tolist()
# extract library names, which we'll use to control the rest of the workflow
# some libraries are split between multiple SRA accessions
ILLUMINA_LIB_NAMES = metadata_illumina["library_name"].unique().tolist()

rule all:
    input: expand("outputs/khmer/{illumina_lib_name}.fq.gz", illumina_lib_name = ILLUMINA_LIB_NAMES)

rule download_fastq_files:
    output: "inputs/raw/{run_accession}.fq.gz"
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
    params: liblayout = lambda wildcards: metadata_illumina.loc[wildcards.illumina_lib_name, "library_layout"]
    shell:'''
    if [ "{params.liblayout}" == "PAIRED" ]; then
        fastp -i {input} --trim_poly_x --qualified_quality_phred 2 --json {output.json} --html {output.html} --report_title {wildcards.illumina_lib_name} --interleaved_in --stdout | gzip > {output.fq}
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        fastp -i {input} --trim_poly_x --qualified_quality_phred 2 --json {output.json} --html {output.html} --report_title {wildcards.illumina_lib_name} --stdout | gzip > {output.fq}
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
    trim-low-abund.py -V -k 20 -Z 18 -C 2 {input} -o {output} -M 4e9 --diginorm --diginorm-coverage=20 --gzip {input}
    '''
# kmer trimming will expect 
