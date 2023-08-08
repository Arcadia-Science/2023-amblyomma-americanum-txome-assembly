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
RUNACCS = metadata_all["run_accession"].unique().tolist()
# extract library names, which we'll use to control the rest of the workflow
# some libraries are split between multiple SRA accessions
ILLLIBNAMES = metadata_illumina["library_name"].unique().tolist()

rule all:
    input: expand("outputs/raw_combined/{illlibname}.fq.gz", illlibname = ILLLIBNAMES)

rule download_fastq_files:
    output: "inputs/raw/{runacc}.fq.gz"
    conda: "envs/sratools.yml"
    params: outdir = "inputs/raw/"
    shell:'''
    fasterq-dump --split-spot -Z {wildcards.runacc} | gzip > {output}
    ''' 


######################################
## Process & assemble illumina files
######################################

rule combine_by_library_name:
    """
    Some of the input sequences have multiple run accessions (SRR*) associated with a single library.
    This rule combines those run accessions (runacc) into one file by illlibname (illumina library name).
    Since it uses the metadata_illumina file to do this, as well as expanding over the runacc wild card in the input,
    the output wildcard illlibname will only contain illumina library names.
    """
    input: expand("inputs/raw/{runacc}.fq.gz", runacc = RUNACCS)
    output: expand("outputs/raw_combined/{illlibname}.fq.gz", illlibname = ILLLIBNAMES)
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
        for library_name, runaccs in libdict.items():
            library_paths = []
            for runacc in runaccs:
                path = params.indir + runacc + ".fq.gz"
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
    input: "outputs/raw_combined/{illlibname}.fq.gz"
    output:
        json = "outputs/fastp/{illlibname}.json",
        html = "outputs/fastp/{illlibname}.html",
        fq = "outputs/fastp/{illlibname}.fq.gz"
    conda: "envs/fastp.yml"
    params: liblayout = lambda wildcards: metadata_all.loc[wildcards.illlibname, "library_layout"]
    shell:'''
    if [ "{params.liblayout}" == "PAIRED" ]; then
        fastp -i {input} --trim_poly_x --qualified_quality_phred 2 --json {output.json} --html {output.html} --report_title {wildcards.illlibname} --interleaved_in --detect_adapter_for_pe --stdout | gzip > {output.fq}
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        fastp -i {input} --trim_poly_x --qualified_quality_phred 2 --json {output.json} --html {output.html} --report_title {wildcards.illlibname} --stdout | gzip > {output.fq}
    fi
    '''

rule khmer_kmer_trim_and_normalization:
    """
    K-mer trim and diginorm (digital normalization, or just normalization) according to the eelpond protocol/elvers.
    The oyster river protocol also supports removal of erroneous k-mers through similar methods.
    """
    input: "outputs/fastp/{illlibname}.fq.gz"
    output: "outputs/khmer/{illlibname}.fq.gz"
    conda: "envs/khmer.yml"
    params: liblayout = lambda wildcards: metadata_all.loc[wildcards.illlibname, "library_layout"]
    shell:'''
    if [ "{params.liblayout}" == "PAIRED" ]; then
        trim-low-abund.py -V -k 20 -Z 18 -C 2 {input} -o - -M 4e9 --diginorm --diginorm-coverage=20 | extract-paired-reads.py --gzip -p {output} # note does not save orphaned pairs
    elif [ "{params.liblayout}" == "SINGLE" ]; then
        trim-low-abund.py -V -k 20 -Z 18 -C 2 {input} -o {output} -M 4e9 --diginorm --diginorm-coverage=20
    fi
    '''
# kmer trimming will expect 
