# Use the fixed config
snakemake -f \
    --cores 5 \
    --use-singularity \
    --singularity-prefix ./tmp/ \
    --singularity-args "--bind $(pwd)" \
    --configfile config_apptainer.yaml \
    --printshellcmds \
    --directory $(pwd) \
    --rerun-incomplete