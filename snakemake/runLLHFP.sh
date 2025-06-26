#!/bin/bash

source ~/.bash_profile  # initialize conda/mamba based on your system (can be ~/.bashrc, ~/.profile, or just conda init)
conda activate snakemake8


#screen -L -S flapro ./snakemake8_sge.sh config.yaml 10 #uncomment if running on sge cluster

screen -L -S flapro ./snakemake8_local.sh config.yaml 5 #uncomment if running locally
