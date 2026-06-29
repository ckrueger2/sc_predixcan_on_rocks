#!/bin/bash

#command
usage() {
    echo "Usage: $0 --oid <OID> --cell <CELL> --out_dir <OUT_DIR>" 
    exit 1
}

#command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --oid)
            OID=$2
            shift 2
            ;;
        --cell)
            CELL=$2
            shift 2
            ;;
        --out_dir)
            OUT_DIR=$2
            shift 2
            ;;
        *)
            echo "unknown flag: $1"
            u
sage
            ;;
    esac
done

#check for required arguments
if [[ -z "$OID" || -z "$CELL" || -z "$OUT_DIR" ]]; then
    echo "ERROR: Missing required arguments"
    usage
    exit 1
fi

#script path
REPO=/home/ckrueger2/sc_pwas_in_ukb

#set up S-PrediXcan environment
#bash "$REPO/00install_predixcan.sh"

#activate conda
source $HOME/miniconda3/bin/activate

#create environment with compatible versions (version numbers may need to be changed with future updates)
if ! conda env list | grep -qw imlabtools; then
    #if it doesn't exist, create it
    conda create -n imlabtools python=3.8 numpy pandas scipy -y
fi

#activate imlabtools
if conda activate imlabtools; then
    echo "Successfully activated imlabtools environment"
fi

output_file="${OUT_DIR}/ukb_${OID}_${CELL}_predixcan_output.csv"

if [ -f "$output_file" ]; then
    echo "WARNING: Output file $output_file already exists. Replacing..."
    rm -f "$output_file"
fi

     
echo "Running PWAS for $OID"

#run s-predixcan - continue even if it fails
export PYTHONNOUSERSITE=1
if python $REPO/01run_predixcan.py --oid "$OID" --cell "$CELL" --out_dir "$OUT_DIR"; then
    echo ""

else
    echo "ERROR: S-PrediXcan failed for $OID (exit code $?)"
fi

#deactivate imlabtools
conda deactivate

#how to view generated PNG files
echo "DONE: View S-PrediXcan file in ${OUT_DIR}"
