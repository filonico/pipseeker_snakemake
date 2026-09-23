# Snakemake make workflow to run PIPseeker

In this directory you can find a Snakemake workflow to run the [PIPseeker™](https://www.fluentbio.com/products/pipseeker-software-for-data-analysis/) software and map single-cell RNA-seq reads derived from the PIPseq chemistry. PIPseq is an emulsion-based single cell encapsulation and barcoding protocol, thus microfluidics free.

The PIPseq chemistry was originally developed by [Fluent BioSciences](https://www.fluentbio.com/), and PIPseeker™ is the freely-available software used to process single-cell RNA-seq data produced with such chemistry.

The methodology is described in this publication:

> Clark et al. 2023. **Microfluidics-free single-cell genomics with templated emulsification**. *Nature biotechnology*, *41*(11), 1557-1566.

The PIPseq chemistry is now held by [Illumina](https://www.illumina.com/techniques/sequencing/rna-sequencing/ultra-low-input-single-cell-rna-seq/pip-seq-chemistry.html), which developed another software to process and map reads, [DRAGEN Single Cell RNA](https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/dragen-single-cell-rna.html). However, this software is apparently available under a subscription fee only...

> [!WARNING]
> PIPseeker last release was version 3.3. **PIPseeker is not actively maintained anymore**. Though, it remains free, as any other existing bioinformatic software...