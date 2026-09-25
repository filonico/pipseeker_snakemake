#!/bin/bash

snakemake --workflow-profile slurm_profile_config.yaml --configfile config.yaml -s Snakefile -n

snakemake --dag --configfile config.yaml -s Snakefile | dot -Tsvg > dag.svg
