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
    input: expand("outputs/raw_combined/{illlibname}.fq.gz", illlibname = ILLLIBNAMES

rule download_fastq_files:
    output: "inputs/raw/{runacc}.fq.gz"
    conda: "envs/sratoolkit.yml"
    params: outdir = "inputs/raw/"
    shell:'''
    fasterq-dump --split-spot -Z {wildcards.runacc} | gzip > {output}
    ''' 

#    params: liblayout = lambda wildcards: metadata_all.loc[wildcards.runacc, "library_layout"]

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
    output: "outputs/raw_combined/{illlibname}.fq.gz"
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


# with single end & paired end files, 
# kmer trimming will expect 
