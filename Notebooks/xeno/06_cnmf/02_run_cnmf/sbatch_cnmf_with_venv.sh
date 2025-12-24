#!/bin/zsh

#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH -J cNMF
. /etc/profile.d/huji-lmod.sh
bash
source /sci/labs/yotamd/lab_share/avishai.wizel/python_envs/Virtual_env/cnmf_dev/bin/activate

cd "/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/EGFR"

python "./Notebooks/xeno/04_cnmf/02_run_cnmf/cnmf_and_harmony_noTPM_100iter_densify.py" \
  "./Reports/xeno/04_cnmf/02_run_cnmf_1.5"\
  "./Reports/xeno/04_cnmf/01_create_data_for_cnmf/xeno_counts_filtered.h5ad"\



