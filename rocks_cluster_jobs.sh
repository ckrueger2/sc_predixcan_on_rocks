#!/bin/bash
#$ -cwd
#$ -V
#$ -N hep_perip
#$ -o /home/ckrueger2/logs/hepatocyte_periportal/sge_$TASK_ID.out
#$ -j y
#$ -l mem_free=19G
#$ -t 1-2158
#$ -tc 50
#$ -pe smp 5

#assign variables
CELL="hepatocyte_periportal"

OUT_DIR="$HOME/sc_pwas_in_ukb/output/${CELL}_cell"
LOG_DIR="$HOME/logs/${CELL}"

mkdir -p "$LOG_DIR"
mkdir -p "$OUT_DIR"

OID=$(sed -n "${SGE_TASK_ID}p" $HOME/sc_pwas_in_ukb/oids_${CELL}_oom.txt)

#rename logs to include OID
mv "$LOG_DIR/sge_${SGE_TASK_ID}.out" "$LOG_DIR/${OID}_${CELL}.out"

#run
bash "$HOME/sc_pwas_in_ukb/00pwas_wrapper.sh" \
    --oid "$OID" --cell "$CELL" --out_dir "$OUT_DIR"
