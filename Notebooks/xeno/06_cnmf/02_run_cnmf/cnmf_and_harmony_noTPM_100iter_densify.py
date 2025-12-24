# to run, open interactive job and:
# cd to this script path
# python cnmf_and_harmony_noTPM_100iter_densify.py
import sys
import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
from cnmf import cNMF, Preprocess
import seaborn as sns
from scipy.sparse import csr_matrix, csc_matrix
matplotlib.use('Agg')


output_folder = sys.argv[1]
batch_correction_folder= output_folder + '/xeno_batchCorrect/'
xeno_filtered_path = sys.argv[2]
num_of_cores = sys.argv[3]

os.makedirs(output_folder, exist_ok=True)
os.makedirs(batch_correction_folder, exist_ok=True)
print(batch_correction_folder)

adata = sc.read(xeno_filtered_path)
adata.X = csr_matrix(adata.X) #convert to sparse to avoid errors

p = Preprocess(random_seed=14)
print ("harmony...")
(adata_c, adata_tp10k, hvgs) = p.preprocess_for_cnmf(adata, harmony_vars='orig.ident',theta = 2, n_top_rna_genes = 2000,max_scaled_thresh = 50,
normalize_librarysize = True, makeplots=False,
save_output_base= batch_correction_folder + 'xeno') # add normalize_librarysize arg to preprocess.py and propogate to normalize_batchcorrect

print ("cnmf prepering...")                                          
cnmf_obj_corrected = cNMF(output_dir=output_folder, name='BatchCorrected_cnmf')
cnmf_obj_corrected.prepare(counts_fn= batch_correction_folder +'xeno.Corrected.HVG.Varnorm.h5ad', genes_file=batch_correction_folder + '/xeno.Corrected.HVGs.txt',components=np.arange(3,11), seed=14, num_highvar_genes=2000,densify=True)

print ("factorizing...")
if (num_of_cores == '1'):
  cnmf_obj_corrected.factorize(worker_i=0, total_workers=1) # without multiprocessing
else:
  cnmf_obj_corrected.factorize_multi_process(total_workers=int(num_of_cores))
  
print ("combining...")
cnmf_obj_corrected.combine()
print ("done...")
cnmf_obj_corrected.k_selection_plot()
import pickle
f = open(output_folder +'/models_2Kvargenes_corrected_noTPM_cnmf_obj.pckl', 'wb')
pickle.dump(cnmf_obj_corrected, f)
f.close()
