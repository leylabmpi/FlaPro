#!/bin/bash

# user input
if [ "$#" -lt 2 ]; then
    echo "snakemake8_local.sh config.yaml cores ..."
    echo " config.yaml : snakemake config"
    echo " cores : number of cores to use locally"
    echo " ... : additional arguments passed to snakemake"
    exit
fi

# check for snakemake
command -v snakemake >/dev/null 2>&1 || { echo "snakemake is not in your PATH"; exit 1; }

# set args
CONFIG=$1
CORES=$2

# snakemake call for local execution
WORKDIR=`pwd`
snakemake -f \
	  --cores $CORES \
	  --use-conda \
	  --configfile $CONFIG \
	  --printshellcmds \
	  --directory $WORKDIR \
	  "${@:3}"