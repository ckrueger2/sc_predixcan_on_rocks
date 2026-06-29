#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import gzip
import pandas as pd
import glob

def set_args():
    parser = argparse.ArgumentParser(description="run s-predixcan with single cell protein prediction")
    parser.add_argument("--oid", help="OlinkID", required=True)
    parser.add_argument("--cell", help="b, hepatocyte_interzonal, hepatocyte_pericentral, cd4_t", required=True)
    parser.add_argument("--out_dir", help="output directory location", required=True)
    return parser
    
def main():
    parser = set_args()
    args = parser.parse_args(sys.argv[1:])
    
    #define paths
    output = f"{args.out_dir}/ukb_{args.oid}_{args.cell}_predixcan_output.csv"
    
    #python and metaxcan paths
    python_path = sys.executable
    metaxcan_dir = "~/MetaXcan"
    
    #retrieve file 
    matches = glob.glob(f"/home/ckrueger2/euro_amer/*_{args.oid}_*.gz")
    if not matches:
        print(f"ERROR: No file found for {args.oid}")
        sys.exit(1)
    gwas_file = matches[0]
    file_name = os.path.basename(gwas_file)
    
    with gzip.open(gwas_file, 'rt') as f:
      df = pd.read_csv(f, sep="\t")
      df["SNP"] = df["CHROM"].astype(str) + "_" + df["GENPOS"].astype(str) + "_" + df["ALLELE0"] + "_" + df["ALLELE1"] + "_b38"
      tmp_file = gwas_file.replace(".gz", "_tmp.txt")
      df.to_csv(tmp_file, sep="\t", index=False)

    #assign database paths (using DB population)
    model_db_path = f"/home/ckrueger2/sc_pwas_in_ukb/models/{args.cell}/predict_db_Model_training_filtered.db"
    covariance_path = f"/home/ckrueger2/sc_pwas_in_ukb/models/{args.cell}/predict_db_Model_training_filtered.txt.gz"
    
    #verify database files exist
    if not os.path.exists(model_db_path):
        print(f"ERROR: Model database not found: {model_db_path}")
        sys.exit(1)
    if not os.path.exists(covariance_path):
        print(f"ERROR: Covariance file not found: {covariance_path}")
        sys.exit(1)
    
    print(f"\nRunning S-PrediXcan:")
    print(f"  OlinkID:        {args.oid}")
    print(f"  Cell Type:      {args.cell}")
    
    #command without optional parameters
    cmd = f"{python_path} {metaxcan_dir}/software/SPrediXcan.py \
    --gwas_file {tmp_file} \
    --snp_column SNP \
    --effect_allele_column ALLELE1 \
    --non_effect_allele_column ALLELE0 \
    --beta_column BETA \
    --se_column SE \
    --model_db_path {model_db_path} \
    --covariance {covariance_path} \
    --keep_non_rsid \
    --model_db_snp_key rsid \
    --throw \
    --output_file {output}"

    #execute the S-PrediXcan command
    print("Executing S-PrediXcan...")
    exit_code = os.system(cmd)

    if exit_code != 0:
        print(f"ERROR: SPrediXcan.py failed with exit code {exit_code}")
        sys.exit(exit_code)

    print(f"\nS-PrediXcan analysis completed successfully")
    print(f"Output file: {output}")
    
    if os.path.exists(tmp_file):
        os.remove(tmp_file)
    
if __name__ == "__main__":
    main()
