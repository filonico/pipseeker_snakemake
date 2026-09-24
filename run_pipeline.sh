#!/bin/bash

snakemake --rerun-incomplete --conda-prefix /alice-home/2/f/fn76/.snakemake/conda --workflow-profile slurm_profile_config.yaml --configfile config.yaml -s Snakefile
