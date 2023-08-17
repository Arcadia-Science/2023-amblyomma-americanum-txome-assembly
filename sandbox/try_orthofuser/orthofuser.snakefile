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
    output: "outputs/orthofuser/transrate_full/assemblies.csv"
    singularity: "docker://pgcbioinfo/transrate:1.0.3"
    params:
        outdir="outputs/orthofuser/transrate_full/"
    threads: 6
    shell:'''
    transrate -o {params.outdir} -t {threads} -a {input.assembly} --left {input.r1} --right {input.r2}
    '''

rule separate_contig_names_into_one_orthogroup_per_file:
    """
    1. Find the Orthogroups.txt file and count its lines:
       wc -l: Counts the number of lines in the file.
       The result (number of lines) is stored in the END environment variable.

    2. Generate a sequence of numbers from 1 to END:
       The numbers 1 through END are generated, each separated by a space.
       tr ' ' '\n': This translates spaces to newlines, so each number is on a new line.
       The result is written to a file named list.

    3. Process each line of the list file in parallel:
       For each number in the list file, extract the corresponding line from the Orthogroups.txt file.
       Split that line into multiple lines (one word or sequence of characters per line).
       Remove the first of those lines.
       Write the result to a .groups file in the outputs directory, with the filename being the number from the list file followed by .groups.
    """
    input: 
        orthogroups = "outputs/orthofuser/orthofinder/Orthogroups/Orthogroups.txt"
    output:
    params: outdir = "outputs/orthofuser/tmp"
    shell:'''
    # export the length of the orthogroups file as a shell var
    export END=$(wc -l {input.orthogroups} | awk '{{print $1}}') && \
    # make a file called "list" that has numbers 1:END, one per line
    echo $(eval echo "{{1..$END}}") | tr ' ' '\n' > list
    # # separate the groups out into different files
    cat list | parallel  -j 1 -k "sed -n ''{{}}'p' {input.orthogroups} | tr ' ' '\n' | sed '1d' > {params.outdir}/{{1}}.groups"
    '''

rule add_transrate_info_to_each_orthogroup:
    """
    Using the *.groups files with orthogroup names, add transrate transcript scores from contigs.csv to each.
    New files are output to *.orthout
    """
    input: 
        contigs = "outputs/orthofuser/transrate_full/merged/contigs.csv",
        groups = ""
    shell:'''
    # add the transrate results to the orthofinder groups
    ls -rtd test/groups/* | parallel -j {threads} "grep -wf {{}} {input} > {{1}}.orthout"
    '''

rule select_highest_score_transcript_from_each_orthogroup:
    """
    2. Process .orthout files and save results to good.list:
       This finds all .orthout files, processes each line in these files with awk to find the maximum value in the 14th column, and then appends the corresponding value from the 1st column to good.list. 
    """
    shell:'''
    ls -rtd test/groups/*orthout | parallel -j 1 "awk -F, -v max=0 '{if(\$9>max){want=\$1; max=\$9}}END{print want}'" >> test/good.list
    '''

rule filter_by_name:
    """
    keep only contigs that ended up in good.list
    """

rule cdhitest:
    """
    collapse at 0.98 identity, which should be removing nearly identical transcripts
    """

rule run_diamond_to_rescue_real_genes:
    """
    run diamond on the raw, unprocessed transcriptomes
    """

rule determine_what_is_annotated_with_diamond_but_isnt_in_orthfuse_output:
