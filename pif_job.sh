#!/bin/bash
#SBATCH --job-name=pif_ped
#SBATCH --account=eccc_oth
#SBATCH --array=1-36                 # one task per species (1-indexed line numbers in aou_list.txt)
#SBATCH --cpus-per-task=4            # matches chains=4 / parallel_chains=4 in the R script
#SBATCH --mem-per-cpu=8G                  # adjust based on largest species dataset
#SBATCH --time=6:00:00              # adjust based on expected runtime per species
#SBATCH --partition=standard
#SBATCH --output=/home/acs001/PIF/logs/pif_%A_%a.out
#SBATCH --error=/home/acs001/PIF/logs/pif_%A_%a.err

# --- setup ---
set -euo pipefail

# make sure log dir exists (harmless if it already does)
mkdir -p /home/acs001/PIF/logs
# make sure storage dir exists (harmless if it already does)
mkdir -p /gpfs/fs7/eccc/esrp/cws/acs001/PIF
# make sure output dir exists (harmless if it already does)
mkdir -p /gpfs/fs7/eccc/esrp/cws/acs001/PIF/output
# make sure convergence dir exists (harmless if it already does)
mkdir -p /gpfs/fs7/eccc/esrp/cws/acs001/PIF/convergence


#controls where R sessions create temporary files to ensure large data species can run
#without trying to use the home directory (which has very little memory)
# make sure the temporary directory exists
mkdir -p /gpfs/fs7/eccc/esrp/cws/acs001/PIF/output/temp
# create local .Renviron file to control tmpdir etc.
cat >> ~/.Renviron << 'EOF'
TMPDIR=/gpfs/fs7/eccc/esrp/cws/acs001/PIF/output/temp
TMP=/gpfs/fs7/eccc/esrp/cws/acs001/PIF/output/temp
TEMP=/gpfs/fs7/eccc/esrp/cws/acs001/PIF/output/temp
EOF


# container image and path to the R script it will run
CONTAINER="/home/acs001/PIF/r-bbsbayes2.sif"
RSCRIPT_PATH="/home/acs001/PIF/4_gpsc_fit.R"

# path to the file listing all aou codes, one per line
AOU_LIST="/home/acs001/PIF/aou_list.txt"

# pull the aou code corresponding to this array task
aou=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${AOU_LIST}")

if [ -z "${aou}" ]; then
  echo "ERROR: no aou code found for array task ${SLURM_ARRAY_TASK_ID} in ${AOU_LIST}"
  exit 1
fi

export aou

echo "Running species aou=${aou} on array task ${SLURM_ARRAY_TASK_ID} with ${SLURM_CPUS_PER_TASK} cpus"

# --- run the R script inside the apptainer container ---
# --bind mounts /home/acs001 (and thus BBS_Trends_CWS) into the container at the same path
# --env aou passes the SLURM-exported variable through to the container's environment
apptainer exec \
--bind /home/acs001:/home/acs001,/gpfs/fs7/eccc/esrp/cws/acs001:/gpfs/fs7/eccc/esrp/cws/acs001 \
--env aou="${aou}" \
  "${CONTAINER}" \
  Rscript "${RSCRIPT_PATH}"
