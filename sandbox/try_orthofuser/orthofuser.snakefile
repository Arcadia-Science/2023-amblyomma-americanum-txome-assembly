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
ASSEMBLERS = ["trinity", "rnaspades"]

READS = ['R1', 'R2']

rule all:
    input: "outputs/orthofuser/orthofuser_final.fa"

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

rule filter_by_length:
    """
    Orthofuser filters to nucleotides greater than 200bp
    We'll use 75, since we're using 25 amino acids as our cut off
    """
    input: "outputs/assembly/renamed/{assembly_group}_{assembler}_renamed.fa"
    output: "outputs/assembly/filtered/{assembly_group}_{assembler}_filtered.fa"
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
        illumina = expand("outputs/assembly/filtered/{assembly_group}_{assembler}_filtered.fa", assembly_group = ASSEMBLY_GROUPS, assembler = ASSEMBLERS),
        isoseq = expand("outputs/assembly/isoseq/{isoseq_lib_name}.fa", isoseq_lib_name = ISOSEQ_LIB_NAMES)
    output: "outputs/orthofuser/orthofinder/Orthogroups/Orthogroups.txt"
    params: 
        indir="outputs/assembly/filtered/",
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

rule merge_txomes:
    input: 
        expand("outputs/assembly/filtered/{assembly_group}_{assembler}_filtered.fa", assembly_group = ASSEMBLY_GROUPS, assembler = ASSEMBLERS),
        expand("outputs/assembly/isoseq/{isoseq_lib_name}.fa", isoseq_lib_name = ISOSEQ_LIB_NAMES)
    output: "outputs/assembly/merged/merged.fa"
    shell:'''
    cat {input} > {output}
    '''

rule combine_reads:
    input: expand("../../outputs/fastp_separated_reads/{illumina_lib_name}_{{read}}.fq.gz", illumina_lib_name = ILLUMINA_LIB_NAMES)
    output: "outputs/fastp_separated_reads_combined/{read}.fq.gz"
    shell:'''
    cat {input} > {output}
    '''

rule transrate:
    '''
    TransRate is a tool for reference-free quality assessment of de novo transcriptome assemblies.
    It uses evidence like read mapping rate and length to score each contig.
    orthofuser uses a modified version of transrate,
    but it only uses an updated version of salmon and other bugs, which shouldn't change the functionality
    https://github.com/macmanes-lab/Oyster_River_Protocol/issues/46
    '''
    input: 
        assembly="outputs/assembly/merged/merged.fa",
        reads=expand("outputs/fastp_separated_reads_combined/{read}.fq.gz", read = READS)
    output: "outputs/orthofuser/transrate_full/merged/contigs.csv"
    singularity: "docker://pgcbioinfo/transrate:1.0.3"
    params: outdir="outputs/orthofuser/transrate_full/"
    threads: 28
    shell:'''
    transrate -o {params.outdir} -t {threads} -a {input.assembly} --left {input.reads[0]} --right {input.reads[1]}
    '''

rule get_contig_name_w_highest_transrate_score_for_each_orthogroup:
    """
    This rule replaces a lot of the bash/awk/grep/thousands of file writing steps in orthofuser/ORP with a python script.
    It writes the name of the highest scoring transcript (contig name) from transrate for each orthogroup.
    """
    input: 
        orthogroups = "outputs/orthofuser/orthofinder/Orthogroups/Orthogroups.txt",
        transrate = "outputs/orthofuser/transrate_full/merged/contigs.csv"
    output: "outputs/orthofuser/orthomerged/good.list"
    shell:'''
    python scripts/get_contig_name_w_highest_transrate_score_for_each_orthogroup.py {input.orthogroups} {input.transrate} {output}    
    '''

rule filter_by_name:
    """
    keep only contigs that ended up in good.list
    """
    input:
        fa="outputs/assembly/merged/merged.fa",
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
        fa = "outputs/assembly/renamed/{assembly_group}_{assembler}_renamed.fa",
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
        fa="outputs/assembly/merged/merged.fa",
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
    threads: 4
    shell:'''
    cd-hit-est -M 5000 -T {threads} -c .98 -i {input} -o {output}
    '''
