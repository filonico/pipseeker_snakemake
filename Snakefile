# before starting, make sure to extract 00_input/mnemiopsis_data.tar.gz and 00_input/ppil_genome.tar.gz

import glob

configfile: "config.yaml"
# profile directory: ~/.config/snakemake/profiles/

valid_chemistry = ["v3", "v4", "V"]

#####################
#     FUNCTIONS     #
#####################

def check_config():
    if not isinstance(config.get("libraries"), dict) or not config["libraries"]:
        raise WorkflowError("Please specify a non-empty \"library\" mapping in the config.yaml file.")
    if not isinstance(config.get("pipseq_chemistry"), str) and config["pipseq_chemistry"] not in valid_chemistry:
        raise WorkflowError("Please specify a valid \"chemistry\" in the config.yaml file. Choose among \"v3\", \"v4\", and \"V\".")
    if not isinstance(config.get("genome"), dict) or not config["genome"]:
        raise WorkflowError("Please specify a non-empty \"genome\" mapping in the config.yaml file.")
    if not isinstance(config.get("softwares"), dict) or not config["softwares"]:
        raise WorkflowError("Please specify a non-empty \"software\" mapping in the config.yaml file.")
    
def get_R1(wildcards):
    files = glob.glob(config["libraries"][wildcards.library] + "/*_R1_001.fastq.gz")
    return files[0]

def get_R2(wildcards):
    files = glob.glob(config["libraries"][wildcards.library] + "/*_R2_001.fastq.gz")
    return files[0]


##########################
#     INITIALIZATION     #
##########################

check_config()

genome_ref_outdir = config["genome_ref"]
trimmed_reads_outdir = config["trimmed_reads"]
pipseeker_1st_outdir = config["pipseeker_1st"]
geneExt_output_outdire = config["geneExt_output"]
pipseeker_2nd_outdir = config["pipseeker_2nd"]

log_dir = config["log_dir"]


####################
#     RULE ALL     #
####################

rule all:
    input:
        expand(pipseeker_1st_outdir + "/{library}/star_out.bam", library = config["libraries"])


####################################
#     PREPARE GENOME REFERENCE     #
####################################

rule standardize_gtf:
    input:
        gtf = config["genome"]["gtf"]
    output:
        gff = temp(genome_ref_outdir + "/standardized_agat.gff"),
        gtf = genome_ref_outdir + "/standardized_agat.gtf"
    log:
        stdout = log_dir + "/01_agat.stdout",
        stderr = log_dir + "/01_agat.stderr"
    conda:
        "envs/agat_env.yaml"
    shell:
        """
        agat_convert_sp_gxf2gxf.pl -g {input.gtf} -o {output.gff} &&
        agat_convert_sp_gff2gtf.pl --gff {output.gff} -o {output.gtf} \
            > {log.stdout} 2> {log.stderr}
        """

rule build_reference_file:
    input:
        genome = config["genome"]["fasta"],
        annotation = genome_ref_outdir + "/standardized_agat.gtf"
    output:
        folder = directory(genome_ref_outdir + "/indexed_genome")
    log:
        stdout = log_dir + "/02_build_ref.stdout",
        stderr = log_dir + "/02_build_ref.stderr"
    params:
        pipseeker = config["softwares"]["pipseeker"]
    shell:
        """
        {params.pipseeker} buildmapref \
            --fasta {input.genome} \
            --gtf {input.annotation} \
            --output-path {output.folder} \
            > {log.stdout} 2> {log.stderr}
        """

######################
#     TRIM READS     #
######################

rule trim_reads:
    input:
        R1 = get_R1,
        R2 = get_R2
    output:
        report_R1 = trimmed_reads_outdir + "/{library}/{library}_R1_trimReport.json",
        trimmed_R1 = trimmed_reads_outdir + "/{library}/{library}_R1_001_trimmed.fastq.gz",
        report_R2 = trimmed_reads_outdir + "/{library}/{library}_R2_trimReport.json",
        trimmed_R2 = trimmed_reads_outdir + "/{library}/{library}_R2_001_trimmed.fastq.gz"
    log:
        stdout = log_dir + "/03_{library}_cutadapt.stdout",
        stderr = log_dir + "/03_{library}_cutadapt.stderr"
    conda:
        "envs/cutadapt_env.yaml"
    shell:
        """
        cutadapt --json {output.report_R1} -j 0 -u -78 -o {output.trimmed_R1} {input.R1} &&
        cutadapt --json {output.report_R2} -j 0 -u -105 -o {output.trimmed_R2} {input.R2} \
        > {log.stdout} 2> {log.stderr}
        """


#####################
#     MAP READS     #
#####################

# snakemake --until map_reads --profile ~/.config/snakemake/profiles/slurm_profile.yaml --jobs 3 --rerun-triggers mtime
rule map_reads:
    input:
        R1 = trimmed_reads_outdir + "/{library}/{library}_R1_001_trimmed.fastq.gz",
        R2 = trimmed_reads_outdir + "/{library}/{library}_R2_001_trimmed.fastq.gz",
        genome_index = genome_ref_outdir + "/indexed_genome"
    output:
        bam_file = pipseeker_1st_outdir + "/{library}/star_out.bam"
    log:
        stdout = log_dir + "/{library}/{library}_pipseeker.stdout",
        stderr = log_dir + "/{library}/{library}_pipseeker.stderr"
    params:
        directory = pipseeker_1st_outdir + "/{library}",
        fastq_dir = trimmed_reads_outdir + "/{library}/{library}",
        chemistry = config["pipseq_chemistry"],
        pipseeker = config["softwares"]["pipseeker"]
    shell:
        """
        {params.pipseeker} full \
            --fastq {params.fastq_dir} \
            --star-index-path {input.genome_index} \
            --chemistry {params.chemistry} \
            --output-path {params.directory} \
            --description {wildcards.library} \
            --retain-barcoded-fastqs \
            --sorted-bam \
            > {log.stdout} 2> {log.stderr}
        """

'''
rule sort_bams:
    input:
        unsorted_bam = pipseeker_1st_outdir + "/{library}/star_out.bam"
    output:
        sorted_bam = pipseeker_1st_outdir + "/{library}/star_out_sorted.bam"
    log:
        stdout = pipseeker_1st_outdir + "/{library}/{library}_samtools_sort.stdout",
        stderr = pipseeker_1st_outdir + "/{library}/{library}_samtools_sort.stderr"
    conda:
        "envs/samtools_env.yaml"
    shell:
        """
        samtools sort {input.unsorted_bam} -o {output.sorted_bam} \
            > {log.stdout} 2> {log.stderr}
        """

rule merge_bams:
    input:
        sorted_bam = expand("03_pipseeker_1st_round/{library}/star_out_sorted.bam", library=config["libraries"])
    output:
        merged_bam = "03_pipseeker_1st_round/merged_bam.bam"
    log:
        stdout = "03_pipseeker_1st_round/samtools_merge.stdout",
        stderr = "03_pipseeker_1st_round/samtools_merge.stderr"
    conda:
        "../envs/samtools_env.yaml"
    shell:
        """
        samtools merge {output.merged_bam} {input.sorted_bam} \
            > {log.stdout} 2> {log.stderr}
        """


#######################
#     RUN GENEEXT     #
#######################

rule run_geneExt:
    input:
        annotation = "01_genome_ref/Ppil_agat.gtf",
        merged_bam = "03_pipseeker_1st_round/merged_bam.bam"
    output:
        annotation_geneExt = "04_genome_ref_geneExt/Ppil_agat_geneExt.gtf"
    log:
        stdout = "04_genome_ref_geneExt/geneExt.stdout",
        stderr = "04_genome_ref_geneExt/geneExt.stderr"
    params:
        geneExt = config["softwares"]["geneExt"]
    conda:
        "../envs/geneExt_env.yaml"
    shell:
        """
        python {params.geneExt} \
            -g {input.annotation} \
            -b {input.merged_bam} \
            -o {output.annotation_geneExt} \
            --force \
            --orphan \
            > {log.stdout} 2> {log.stderr}
        """

###########################
#     MAP READS AGAIN     #
###########################

rule build_reference_file_geneExt:
    input:
        genome = config["genome"]["fasta"],
        annotation = "04_genome_ref_geneExt/Ppil_agat_geneExt.gtf"
    output:
        folder = directory("04_genome_ref_geneExt/indexed_genome")
    log:
        stdout = "04_genome_ref_geneExt/build_ref.stdout",
        stderr = "04_genome_ref_geneExt/build_ref.stderr"
    params:
        pipseeker = config["softwares"]["pipseeker"]
    shell:
        """
        {params.pipseeker} buildmapref \
            --fasta {input.genome} \
            --gtf {input.annotation} \
            --output-path {output.folder} \
            > {log.stdout} 2> {log.stderr}
        """

# snakemake --until map_reads_geneExt --profile ~/.config/snakemake/profiles/slurm_profile.yaml --jobs 3 --rerun-triggers mtime
rule map_reads_geneExt:
    input:
        R1 ="02_trimmed_reads/{library}/{library}_R1_001_trimmed.fastq.gz",
        R2 = "02_trimmed_reads/{library}/{library}_R2_001_trimmed.fastq.gz",
        genome_index = "04_genome_ref_geneExt/indexed_genome"
    output:
        directory = directory("05_pipseeker_2nd_round/{library}")
    log:
        stdout = "05_pipseeker_2nd_round/{library}/{library}_pipseeker.stdout",
        stderr = "05_pipseeker_2nd_round/{library}/{library}_pipseeker.stderr"
    params:
        fastq_dir = "02_trimmed_reads/{library}/{library}",
        pipseeker = config["softwares"]["pipseeker"]
    resources:
        runtime = 480,
        mem_mb = 70000,
        cpus_per_task = 20
    shell:
        """
        {params.pipseeker} full \
            --fastq {params.fastq_dir} \
            --star-index-path {input.genome_index} \
            --chemistry V \
            --output-path {output.directory} \
            --description {wildcards.library} \
            --retain-barcoded-fastqs \
            > {log.stdout} 2> {log.stderr}
        """

###########################################
#     FORCE PIPSEEKER CELLS AND MERGE     #
###########################################

# not run
rule force_pipseeker_cells:
    input:
        previous_run = "05_pipseeker_2nd_round/{library}"
    output:
        done = touch("06_pipseeker_merge/.{library}_force_done")
    log:
        stdout = "05_pipseeker_2nd_round/{library}/{library}_pipseeker_force.stdout",
        stderr = "05_pipseeker_2nd_round/{library}/{library}_pipseeker_force.stderr"
    params:
        pipseeker = config["softwares"]["pipseeker"]
    shell:
        """
        {params.pipseeker} cells \
            --previous {input.previous_run} \
            --force-cells 20000 \
            > {log.stdout} 2> {log.stderr}
        """

# not run
rule merge_pipseeker:
    input:
        forceCell_results = expand("06_pipseeker_merge/.{library}_force_done", library=config["libraries"])
    output:
        merged_output = directory("07_pipseeker_merge_batch")
    log:
        stdout = "06_pipseeker_merge_batch/pipseeker_merge.stdout",
        stderr = "06_pipseeker_merge_batch/pipseeker_merge.stderr"
    params:
        pipseeker = config["softwares"]["pipseeker"],
        previous_runs_joined = ",".join(expand("05_pipseeker_2nd_round/{library}", library = config["libraries"])),
        labels = ",".join(config["libraries"])
    shell:
        """
        {params.pipseeker} merge \
            --previous {params.previous_runs_joined} \
            --cell-calling-mode force_20000,force_20000,force_20000 \
            --output-path {output.merged_output} \
            --sample-labels "{params.labels}" \
            --batch {params.labels} \
            > {log.stdout} 2> {log.stderr}
        """

# this is because apparently pipseeker "force_cell" output is missing a \n character at the end of barcodes.tsv.gz
# so we need to add it
rule update_matrices:
    input:
        force_cell = "06_pipseeker_merge/.{library}_force_done"
    output:
        updated_mappings = directory("08_input_for_seurat/force_20000/{library}")
    log:
        stdout = "08_input_for_seurat/{library}.stdout",
        stderr = "08_input_for_seurat/{library}.stderr"
    shell:
        """
        cp -r 05_pipseeker_2nd_round/{wildcards.library}/filtered_matrix/force_20000 {output.updated_mappings} &&
            gunzip {output.updated_mappings}/barcodes.tsv.gz &&
            echo "" >> {output.updated_mappings}/barcodes.tsv &&
            gzip {output.updated_mappings}/barcodes.tsv \
            > {log.stdout} 2> {log.stderr}
        """
'''