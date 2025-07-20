#!/bin/zsh

#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -J cNMF
#SBATCH --mail-user=avishai.wizel@mail.huji.ac.il
#SBATCH --mail-type=END

. /etc/profile.d/huji-lmod.sh
bash
export PATH="/sci/labs/yotamd/lab_share/avishai.wizel/python_envs/miniconda/bin/:$PATH"
eval "$(conda shell.bash hook)"

conda activate cnmf_1.7
cd "/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/EGFR"

python "./Notebooks/xeno/04_cnmf/02_run_cnmf/cnmf_and_harmony_noTPM_100iter_densify.py" \
  "./Reports/xeno/04_cnmf/02_run_cnmf_1.7"\
  "./Reports/xeno/04_cnmf/01_create_data_for_cnmf/xeno_counts_filtered.h5ad"\



