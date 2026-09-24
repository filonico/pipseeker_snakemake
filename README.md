# Snakemake make workflow to run PIPseeker

In this directory you can find a Snakemake workflow to run the [PIPseeker™](https://www.fluentbio.com/products/pipseeker-software-for-data-analysis/) software and map single-cell RNA-seq reads derived from the PIPseq chemistry. PIPseq is an emulsion-based single cell encapsulation and barcoding protocol, thus microfluidics free.

The PIPseq chemistry was originally developed by [Fluent BioSciences](https://www.fluentbio.com/), and PIPseeker™ is the freely-available software used to process single-cell RNA-seq data produced with such chemistry.

The methodology is described in this publication:

> Clark et al. 2023. **Microfluidics-free single-cell genomics with templated emulsification**. *Nature biotechnology*, *41*(11), 1557-1566. doi: [10.1038/s41587-023-01685-z](https://doi.org/10.1038/s41587-023-01685-z)

The PIPseq chemistry is now held by [Illumina](https://www.illumina.com/techniques/sequencing/rna-sequencing/ultra-low-input-single-cell-rna-seq/pip-seq-chemistry.html), which developed another software to process and map reads, [DRAGEN Single Cell RNA](https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/dragen-single-cell-rna.html). However, this software is apparently available under a subscription fee only...

> [!WARNING]
> PIPseeker last release was version 3.3. **PIPseeker is not actively maintained anymore**. Though, it remains free, as any other existing bioinformatic software...

# How to run this pipeline
1. **Clone this repository** onto your working directory.

```bash
git clone https://github.com/filonico/pipseeker_snakemake.git
```

2. **Download and install PIPseeker v3.3.0**, which does not currently have any conda distribution. It can be obtained from this webpage [fluentbio.com/resources/pipseeker-downloads/](https://www.fluentbio.com/resources/pipseeker-downloads/). You may be asked to insert your details.

3. **Download the GeneExt executable** from [GitHub](https://github.com/sebepedroslab/GeneExt). You don't need to create the corresponding conda environment, this will be handled by snakemake.

4. **Have your PIPseq scRNA-seq raw reads and annotated genome assembly ready**. Read fastq files need to be splitted in R1 and R2 reads. The genome needs to have a least the assembly fasta file and the annotation gtf/gff file.

5. **Edit both the [`config.yaml`](./config.yaml) and the [`slurm_profile_config.yaml`](./slurm_profile_config.yaml) files** with your data.

6. **Execute a Snakemake dry run**.

```bash
snakemake --rerun-incomplete \
    --workflow-profile slurm_profile_config.yaml \
    --configfile config.yaml \
    -s Snakefile
```