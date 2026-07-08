# scPrediXcan on Loyola Medical School ROCKS cluster
### The following scripts can be used to run single-cell prediction models using S-PrediXcan on a ROCKS cluster

#### Must be edited:
- **rocks_cluster_jobs_USER.sh:** Runs 00pwas_wrapper using scheduler. Change file paths. Use rocks_cluster_jobs as an example, if needed
- **rocks_combine.sh:** Runs 02format_output.R using scheduler. Change file paths. ex run. qsub ~/sc_pwas_in_ukb/rocks_combine.sh hepatocyte_interzonal
- **02format_output.R:** Combines all output files, filtering by p<0.05 and bonferroni correction
#### Should not need to be edited (until my account is deleted):
- **00install_predixcan.sh:** Installs miniconda and MetaXcan repo onto the cluster, then builds the imlabtools conda environment
- **00pwas_wrapper.sh:** Wrapper that runs PWAS with scPrediXcan
- **01run_predixcan.py:** Executes scPrediXcan via wrapper
- **oids.txt:** List of all UK Biobank OIDs available
#### Other:
- **missing_oid_outputs.R:** Can be used to check for files present in logs, but not scPrediXcan output or files present in logs, but not full oid list. Change file paths. If used, change OID file in rocks_cluster_jobs_USER.sh to new file name.
- **03mapping.R:** Maps gene coordinates, classifies cis same/different and trans from Wittich et al. 2024, calculates FDR significance (I did this script on compbio1)
#### Running scripts:
- qsub script.sh (to run a job)
- qstat (see jobs running)
- qdel (kill/remove job)
- qstat -j <job_id> | grep -A3 error
- qstat -u USERNAME
#### Copy files to other servers (example whole folder, individual file(s)):
- scp -r /home/ckrueger2/sc_pwas_in_ukb/output/hepatocyte_interzonal_cell/ claudia@compbio1.cs.luc.edu:/home/claudia/sc_output/
- scp /home/ckrueger2/sc_pwas_in_ukb/filtered_output/* claudia@compbio1.cs.luc.edu:/home/claudia/sc_output_filtered
