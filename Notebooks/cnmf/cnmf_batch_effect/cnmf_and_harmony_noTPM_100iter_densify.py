# to run, open interactive job and:
# source /sci/labs/yotamd/lab_share/avishai.wizel/python_envs/Virtual_env/cnmf_dev/bin/activate.csh
# cd "/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/EGFR/Notebooks/cnmf/cnmf_batch_effect/"
# python cnmf_and_harmony_noTPM_100iter_densify.py

import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
import scanpy as sc
from cnmf import cNMF, Preprocess
import seaborn as sns

adata = sc.read("/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/EGFR/Data/cnmf_V2/xeno_counts_filtered.h5ad")
p = Preprocess(random_seed=14)
print ("harmony...")
(adata_c, adata_tp10k, hvgs) = p.preprocess_for_cnmf(adata, harmony_vars='orig.ident',theta = 2, n_top_rna_genes = 2000,max_scaled_thresh = 50, quantile_thresh = None,librarysize_targetsum = None, makeplots=True,save_output_base='/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/EGFR/Data/cnmf_V2/noTPM_100iter/xeno_batchCorrect')

print ("cnmf prepering...")                                          
cnmf_obj_corrected = cNMF(output_dir='/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/EGFR/Data/cnmf_V2/noTPM_100iter/', name='BatchCorrected_cnmf')
cnmf_obj_corrected.prepare(counts_fn='/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/EGFR/Data/cnmf_V2/noTPM_100iter/xeno_batchCorrect.Corrected.HVG.Varnorm.h5ad', genes_file='/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/EGFR/Data/cnmf_V2/noTPM_100iter/xeno_batchCorrect.Corrected.HVGs.txt',components=np.arange(3,11), seed=14, num_highvar_genes=2000,densify=True)

print ("factorizing...")
cnmf_obj_corrected.factorize(worker_i=0, total_workers=1)
print ("combining...")
cnmf_obj_corrected.combine()
print ("done...")
import pickle
f = open('/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/EGFR/Data/cnmf_V2/noTPM_100iter/models_2Kvargenes_corrected_noTPM_cnmf_obj.pckl', 'wb')
pickle.dump(cnmf_obj_corrected, f)
f.close()
