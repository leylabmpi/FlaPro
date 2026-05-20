# Use the fixed config
snakemake -f \
    --cores 5 \
    --use-apptainer \
    --apptainer-prefix ./tmp/ \
    --apptainer-args "--bind $(pwd)" \
    --configfile config_apptainer_orig.yaml \
    --printshellcmds \
    --directory $(pwd) \
    --rerun-incomplete