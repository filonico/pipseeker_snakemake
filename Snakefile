# before starting, make sure to extract 00_input/mnemiopsis_data.tar.gz and 00_input/ppil_genome.tar.gz

import glob

configfile: "config.yaml"
# profile directory: ~/.config/snakemake/profiles/

valid_chemistry = ["v3", "v4", "V"]


##########################################
#     FUNCTIONS TO RUN SAFETY CHECKS     #
##########################################

def check_config():
    if not isinstance(config.get("libraries"), dict) or not config["libraries"]:
        raise WorkflowError("Please specify a non-empty \"library\" mapping in the config.yaml file.")
    if config["pipseq_chemistry"] not in valid_chemistry:
        raise WorkflowError("Please specify a valid \"chemistry\" in the config.yaml file. Choose among \"v3\", \"v4\", and \"V\".")
    if not isinstance(config.get("genome"), dict) or not config["genome"]:
        raise WorkflowError("Please specify a non-empty \"genome\" mapping in the config.yaml file.")
    if not isinstance(config.get("softwares"), dict) or not config["softwares"]:
        raise WorkflowError("Please specify a non-empty \"software\" mapping in the config.yaml file.")
    force_cells = config.get("force_cells", [])
    if force_cells and not isinstance(force_cells, (str, int, list)):
        raise WorkflowError("\"force_cells\" must be a single value or a list of values, if present.")
    
def get_R1(wildcards):
    files = sorted(glob.glob(config["libraries"][wildcards.library] + "/*_R1_001.fastq.gz"))
    if len(files) != 1:
        raise WorkflowError(f"Expected exactly one R1 file for library '{wildcards.library}', found {len(files)}: {files}")
    return files[0]

def get_R2(wildcards):
    files = sorted(glob.glob(config["libraries"][wildcards.library] + "/*_R2_001.fastq.gz"))
    if len(files) != 1:
        raise WorkflowError(f"Expected exactly one R2 file for library '{wildcards.library}', found {len(files)}: {files}")
    return files[0]


##########################
#     INITIALIZATION     #
##########################

check_config()

genome_ref_outdir = config["genome_ref"]
trimmed_reads_outdir = config["trimmed_reads"]
pipseeker_1st_outdir = config["pipseeker_1st"]
geneExt_output_outdir = config["geneExt_output"]
pipseeker_2nd_outdir = config["pipseeker_2nd"]
pipseeker_force_outdir = config["pipseeker_force"]

log_dir = config["log_dir"]

# turn force_cells to a list of strings; if empty, skip force-cells-related steps
_force_cells_raw = config.get("force_cells", [])
if isinstance(_force_cells_raw, (str, int)):
    force_cells_list = [str(_force_cells_raw)]
elif isinstance(_force_cells_raw, list):
    force_cells_list = [str(v) for v in _force_cells_raw]
else:
    force_cells_list = []


####################
#     RULE ALL     #
####################

def get_rullAll_inputs():
    # if force-cell is specified, then use them to list of inputs
    if force_cells_list:
        inputs = expand(
            pipseeker_force_outdir + "/force_{force_cells}/{library}",
            library = config["libraries"],
            force_cells = force_cells_list
        )
    # else get only pipseeker 2nd round of mapping out dir
    else:
        inputs = expand(pipseeker_2nd_outdir + "/{library}", library=config["libraries"])
    
    return inputs

rule all:
    input:
        get_rullAll_inputs()


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
        folder = directory(genome_ref_outdir + "/indexed_genome"),
        done_placeholder = touch(genome_ref_outdir + "/.indexed_genome_done") # creating an actual file is safer than only a directory
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

rule map_reads:
    input:
        R1 = trimmed_reads_outdir + "/{library}/{library}_R1_001_trimmed.fastq.gz",
        R2 = trimmed_reads_outdir + "/{library}/{library}_R2_001_trimmed.fastq.gz",
        genome_index = genome_ref_outdir + "/indexed_genome"
    output:
        bam_file = pipseeker_1st_outdir + "/{library}/star_out.bam"
    log:
        stdout = log_dir + "/04_{library}_pipseeker.stdout",
        stderr = log_dir + "/04_{library}_pipseeker.stderr"
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
            > {log.stdout} 2> {log.stderr}
        """


# pipseeker has a "--sorted_bam" flag to output sorted bam files, but it doesn't work (unrecognised arguemnt)
rule sort_bams:
    input:
        unsorted_bam = pipseeker_1st_outdir + "/{library}/star_out.bam"
    output:
        sorted_bam = pipseeker_1st_outdir + "/{library}/star_out_sorted.bam"
    log:
        stdout = log_dir + "/05_{library}_samtools_sort.stdout",
        stderr = log_dir + "/05_{library}_samtools_sort.stderr"
    conda:
        "envs/samtools_env.yaml"
    shell:
        """
        samtools sort -@ {resources.cpus_per_task} {input.unsorted_bam} -o {output.sorted_bam} \
            > {log.stdout} 2> {log.stderr}
        """

rule merge_bams:
    input:
        sorted_bam = expand(pipseeker_1st_outdir + "/{library}/star_out_sorted.bam", library = config["libraries"])
    output:
        merged_bam = pipseeker_1st_outdir + "/merged_bam.bam"
    log:
        stdout = log_dir + "/06_samtools_merge.stdout",
        stderr = log_dir + "/06_samtools_merge.stderr"
    conda:
        "envs/samtools_env.yaml"
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
        annotation = genome_ref_outdir + "/standardized_agat.gtf",
        merged_bam = pipseeker_1st_outdir + "/merged_bam.bam"
    output:
        annotation_geneExt = geneExt_output_outdir + "/Ppil_agat_geneExt.gtf"
    log:
        stdout = log_dir + "/07_geneExt.stdout",
        stderr = log_dir + "/07_geneExt.stderr"
    params:
        geneExt = config["softwares"]["geneExt"]
    conda:
        "envs/geneExt_env.yaml"
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
        annotation = geneExt_output_outdir + "/Ppil_agat_geneExt.gtf"
    output:
        folder = directory(geneExt_output_outdir + "/indexed_genome")
    log:
        stdout = log_dir + "/08_build_ref_geneExt.stdout",
        stderr = log_dir + "/08_build_ref_geneExt.stderr"
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

rule map_reads_geneExt:
    input:
        R1 = trimmed_reads_outdir + "/{library}/{library}_R1_001_trimmed.fastq.gz",
        R2 = trimmed_reads_outdir + "/{library}/{library}_R2_001_trimmed.fastq.gz",
        genome_index = geneExt_output_outdir + "/indexed_genome"
    output:
        directory = directory(pipseeker_2nd_outdir + "/{library}"),
        done_placeholder = touch(pipseeker_2nd_outdir + "/.{library}_done") # creating an actual file is safer than only a directory
    log:
        stdout = log_dir + "/09_{library}_pipseeker.stdout",
        stderr = log_dir + "/09_{library}_pipseeker.stderr"
    params:
        fastq_dir = trimmed_reads_outdir + "/{library}/{library}",
        pipseeker = config["softwares"]["pipseeker"]
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


#######################################
#     FORCE PIPSEEKER FORCE CELLS     #
#######################################

rule force_pipseeker_cells:
    input:
        previous_run = pipseeker_2nd_outdir + "/{library}"
    output:
        # pipseeker cells adds files to the previously computed run in place,
        # so this placeholder tells snakemake the step has run for this {force_cells} value
        done_placeholder = touch(pipseeker_force_outdir + "/.{library}_force_{force_cells}_done") 
    log:
        stdout = log_dir + "/10_{library}_force_{force_cells}_pipseeker_force.stdout",
        stderr = log_dir + "/10_{library}_force_{force_cells}_pipseeker_force.stderr"
    params:
        pipseeker = config["softwares"]["pipseeker"]
    shell:
        """
        {params.pipseeker} cells \
            --previous {input.previous_run} \
            --force-cells {wildcards.force_cells} \
            > {log.stdout} 2> {log.stderr}
        """

# pipseeker's "force_cells" output is missing a trailing \n in barcodes.tsv.gz, so we add it
rule update_matrices:
    input:
        force_done = pipseeker_force_outdir + "/.{library}_force_{force_cells}_done"
    output:
        updated_mappings = directory(pipseeker_force_outdir + "/force_{force_cells}/{library}"),
        done_placeholder = touch(pipseeker_force_outdir + "/force_{force_cells}/.{library}_done") # creating an actual file is safer than only a directory
    log:
        stdout = log_dir + "/11_{library}_force_{force_cells}_pipseeker_update_matrices.stdout",
        stderr = log_dir + "/11_{library}_force_{force_cells}_pipseeker_update_matrices.stderr"
    params:
        previous_run = pipseeker_2nd_outdir + "/{library}"
    shell:
        """
        cp -r {params.previous_run}/filtered_matrix/force_{wildcards.force_cells} {output.updated_mappings} &&
            gunzip {output.updated_mappings}/barcodes.tsv.gz &&
            echo "" >> {output.updated_mappings}/barcodes.tsv &&
            gzip {output.updated_mappings}/barcodes.tsv \
            > {log.stdout} 2> {log.stderr}
        """

'''
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
'''