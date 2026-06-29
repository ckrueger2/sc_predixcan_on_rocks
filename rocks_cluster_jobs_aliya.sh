#!/bin/bash
#$ -cwd
#$ -V
#$ -N job_submission_name
#$ -o /home/ahaas4/logs/sge_$TASK_ID.out
#$ -j y
#$ -l mem_free=19G
#$ -t 1-2896
#$ -tc 50

#note: make sure to change job_submission_name and that -t range is the correct number of proteins in oids.txt file used

#assign variables
CELL="cholangiocyte"

OUT_DIR="/home/ahaas4/spredixcan/output/${CELL}_cell"
LOG_DIR="/home/ahaas4/logs"

mkdir -p "$LOG_DIR"
mkdir -p "$OUT_DIR"

OID=$(sed -n "${SGE_TASK_ID}p" /home/ckrueger2/sc_pwas_in_ukb/oids.txt)

#rename logs to include OID
mv "$LOG_DIR/sge_${SGE_TASK_ID}.out" "$LOG_DIR/${OID}_${CELL}.out"

#run
bash "/home/ckrueger2/sc_pwas_in_ukb/00pwas_wrapper.sh" \
    --oid "$OID" --cell "$CELL" --out_dir "$OUT_DIR"
