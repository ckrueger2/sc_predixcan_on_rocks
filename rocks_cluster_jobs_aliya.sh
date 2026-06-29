#!/bin/bash
#$ -cwd
#$ -V
#$ -N job_submission_name
#$ -o /home/YOUR_USERNAME/logs/sge_$TASK_ID.out
#$ -j y
#$ -l mem_free=19G
#$ -t 1-2896 #number of jobs
#$ -tc 50 #how many jobs to run at maximum (65 total on server, check other user's runs prior to running)
#$ -pe smp 5 #ensures 1 job per server

#note: make sure to change job_submission_name, output username path, and that -t range is the correct number of proteins in oids.txt file used

#assign variables (change accordingly)
CELL="hepatocyte_periportal"
USER = "ckrueger2"

OUT_DIR="/home/${USER}/sc_pwas_in_ukb/output/${CELL}_cell"
LOG_DIR="/home/${USER}/logs/${CELL}"

mkdir -p "$LOG_DIR"
mkdir -p "$OUT_DIR"

OID=$(sed -n "${SGE_TASK_ID}p" /home/ckrueger2/sc_pwas_in_ukb/oids.txt) #change this path for any re-runs

#rename logs to include OID
mv "$LOG_DIR/sge_${SGE_TASK_ID}.out" "$LOG_DIR/${OID}_${CELL}.out"

#run
bash "/home/ckrueger2/sc_pwas_in_ukb/00pwas_wrapper.sh" \
    --oid "$OID" --cell "$CELL" --out_dir "$OUT_DIR"
