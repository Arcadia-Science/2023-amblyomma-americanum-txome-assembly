# start point -- combined assemblies
# actually orthofuser may expect everything to be in separate files, but all in the same folder

ASSEMBLY_GROUPS = ['petx120midgutfemale', 'petx120sgfemale']

rule all:
    input: 
        "outputs/orthofuser/orthofinder/Orthogroups/Orthogroups.txt",
        "outputs/orthofuser/transrate_full/assemblies.csv"

rule rename_contigs:
    """
    we should talk to austin and figure out what his requirements are for transcript FASTA headers for noveltree
    """
    input: "assemblies/{assembly_group}_rnaspades_soft_filtered_transcripts.fa"
    output: "outputs/assemblies/renamed/{assembly_group}_renamed.fa"
    conda: "envs/bbmap.yml"
    shell:'''
    bbrename.sh in={input} out={output} prefix={wildcards.assembly_group} addprefix=t
    '''

rule filter_by_length:
    """
    Orthofuser filters to nucleotides greater than 200bp
    We'll use 75, since we're using 25 amino acids as our cut off
    """
    input: "outputs/assemblies/renamed/{assembly_group}_renamed.fa"
    output: "outputs/assemblies/filtered/{assembly_group}_filtered.fa"
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
    input: expand("outputs/assemblies/filtered/{assembly_group}_filtered.fa", assembly_group = ASSEMBLY_GROUPS)
    output: "outputs/orthofuser/orthofinder/Orthogroups/Orthogroups.txt"
    params: 
        indir="outputs/assemblies/filtered/",
        outdir = "outputs/orthofuser/orthofinder",
        outdirtmp = "outputs/orthofuser/orthofinder_tmp"
    threads: 7
    conda: "envs/orthofinder.yml"
    shell:'''
    orthofinder -d -I 12 -f {params.indir} -o {params.outdirtmp} -og -t {threads} -a {threads} && mv {params.outdirtmp}/Results*/* {params.outdir}
    '''

rule merge_txomes:
    input: expand("outputs/assemblies/filtered/{assembly_group}_filtered.fa", assembly_group = ASSEMBLY_GROUPS)
    output: "outputs/assemblies/merged/merged.fa"
    shell:'''
    cat {input} > {output}
    '''

rule transrate:
    '''
    orthofuser uses a modified version of transrate,
    but it only uses an updated version of salmon and other bugs, which shouldn't change the functionality
    https://github.com/macmanes-lab/Oyster_River_Protocol/issues/46
    '''
    input: 
        assembly="outputs/assemblies/merged/merged.fa",
        r1="left.fq.gz",
        r2="right.fq.gz"
    output: "outputs/orthofuser/transrate_full/merged/contigs.csv"
    singularity: "docker://pgcbioinfo/transrate:1.0.3"
    params:
        outdir="outputs/orthofuser/transrate_full/"
    threads: 6
    shell:'''
    transrate -o {params.outdir} -t {threads} -a {input.assembly} --left {input.r1} --right {input.r2}
    '''

rule get_contig_name_w_highest_transrate_score_for_each_orthogroup:
    """
    This rule replaces a lot of the bash/awk/grep/thousands of file writing steps in orthofuser/ORP with a python script.
    It write the name of the highest scoring transcript (contig name) from transrate for each orthogroup.
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
        fa="outputs/assemblies/merged/merged.fa",
        lst="outputs/orthofuser/orthomerged/good.list"
    output: "outputs/orthofuser/orthomerged/orthomerged.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk subseq {input.fa} {input.lst} > {output}
    '''
    

rule run_diamond_on_orthomerged_txome:
    """
    """
    input: 
        fa ="outputs/orthofuser/orthomerged/orthomerged.fa",
        db = "" # swissprot
    output: "outputs/orthofuser/diamond/orthomerged.diamond.txt"
    conda: "envs/diamond.yml"
    threads: 7
    shell:'''
    diamond blastx --quiet -p {threads} -e 1e-8 --top 0.1 -q {input.fa} -d {input.db} -o {output}
    '''

rule run_diamond_to_rescue_real_genes:
    """
    run diamond on the raw, unprocessed transcriptomes
    """

rule determine_what_is_annotated_with_diamond_but_isnt_in_orthfuse_output:

rule cdhitest:
    """
    collapse at 0.98 identity, which should be removing nearly identical transcripts
    """
    input: "outputs/orthofuser/orthomerged/orthomerged.fa"
    output: "outputs/orthofuser/orthomerged/ORP_intermediate.fa"
    conda: "envs/cd-hit.yml"
    threads: 1
    shell:'''
    cd-hit-est -M 5000 -T {threads} -c .98 -i {input} -o {output}
    '''
