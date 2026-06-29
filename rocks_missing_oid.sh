#!/bin/bash

source ~/.bashrc
conda activate imlabtools

/home/ckrueger2/miniconda3/envs/imlabtools/bin/Rscript /home/ckrueger2/sc_pwas_in_ukb/missing_oid_outputs.R --cell_type "$1"

conda deactivate
