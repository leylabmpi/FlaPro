#!/bin/bash

source ~/.bash_profile
conda activate snakemake8


#screen -L -S flapro ./snakemake8_sge.sh config.yaml 10 #uncomment if running on sge cluster

screen -L -S flapro ./snakemake8_local.sh config.yaml 5 #uncomment if running locally
