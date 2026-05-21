#!/bin/bash
# Requires Snakemake 8. Ensure it is available before running
# (e.g., conda activate snakemake8, module load snakemake/8.x, or equivalent)

snakemake -f \
    --cores 5 \
    --use-apptainer \
    --apptainer-prefix ./tmp/ \
    --apptainer-args "--bind $(pwd)" \
    --configfile config_apptainer.yaml \
    --printshellcmds \
    --directory $(pwd) \
    --rerun-incomplete