library(magrittr)
library(stringr)
library(stringi)

# env path:
cnmf_conda_env_path = "/sci/labs/yotamd/lab_share/avishai.wizel/python_envs/miniconda/envs/cnmf_1.7"

#input data:
patients_count_matrix = "./input_data/patients_raw/fc.txt.gz"
xeno_counts_dir = "./input_data/xeno_raw/"
hif_targets= "./input_data/Pathways/pathways_from_papers/HIF_targets_Lombardi_PMC9869179.txt"
msigDB_hallmarks = "./input_data/Pathways/h.all.v2025.1.Hs.symbols.gmt"
msigdb = "./input_data/Pathways/msigdb.v2025.1.Hs.symbols.RDS"
pipeline = list()
####################################### Preprocess ####################################################

pipeline[["xeno_preprocess"]] = list(
  input = list(
    script = "./Notebooks/xeno/01_preprocess.Rmd"
  ),
  output = list(
    report ="./Reports/xeno/01_preprocess/01_preprocess.html",
    xeno = "./Reports/xeno/01_preprocess/xeno.qs2"),
  params = list(xeno_counts_dir= xeno_counts_dir)
)

pipeline[["patients_preprocess"]] = list(
  input = list(
    script = "./Notebooks/patients/01_preprocess.Rmd",
    count_matrix = patients_count_matrix
  ),
  output = list(
    report ="./Reports/patients/01_preprocess/01_preprocess.html",
    patients = "./Reports/patients/01_preprocess/patients.qs2")
)

####################################### Clustering ####################################################
pipeline[["xeno_clustering"]] = list(
  input = list(
    script = "./Notebooks/xeno/02_clustering.Rmd",
    xeno =  pipeline$xeno_preprocess$output$xeno
  ),
  output = list(
    report ="./Reports/xeno/02_clustering/02_clustering.html"
    )
)


pipeline[["patients_clustering"]] = list(
  input = list(
    script = "./Notebooks/patients/001_clustering.Rmd",
    lung =  pipeline$patients_preprocess$output$patients
  ),
  output = list(
    report ="./Reports/patients/001_clustering/001_clustering.html"
  )
)

####################################### Patients DEG & SIPSIC####################################################


pipeline[["patients_deg"]] = list(
  input = list(
    script = "./Notebooks/patients/01_DEG.Rmd",
    lung =  pipeline$patients_preprocess$output$patients,
    hif_targets = hif_targets,
    genesets = msigDB_hallmarks
  ),
  output = list(
    report ="./Reports/patients/01_DEG/01_DEG.html"
  )
)

pipeline[["patients_run_sipsic"]] = list(
  input = list(
    script ="./Notebooks/patients/02_run_sipsic.Rmd",
    lung =  pipeline$patients_preprocess$output$patients,
    hif_targets= hif_targets,
    genesets= msigDB_hallmarks
  ),
  output = list(
    report ="./Reports/patients/02_run_sipsic/02_run_sipsic.html",
    sipsic_matrix ="./Reports/patients/02_run_sipsic/patients_pathwayScoresMatrix.RDS"
  )
)

pipeline[["patients_sipsic_analysis"]] = list(
  input = list(
    script ="./Notebooks/patients/03_sipsic_analysis.Rmd",
    lung =  pipeline$patients_preprocess$output$patients,
    sipsic_matrix = pipeline$patients_run_sipsic$output$sipsic_matrix
  ),
  output = list(
    report ="./Reports/patients/03_sipsic_analysis/03_sipsic_analysis.html"
  )
)

####################################### Xenogratfs DEG & SIPSIC ####################################################
pipeline[["xenografts_deg"]] = list(
  input = list(
    script = "./Notebooks/xeno/03_DEG.Rmd",
    xeno =  pipeline$xeno_preprocess$output$xeno,
    hif_targets= hif_targets,
    genesets = msigDB_hallmarks
  ),
  output = list(
    report ="./Reports/xeno/03_DEG/03_DEG.html"
  )
)

pipeline[["xenografts_run_sipsic"]] = list(
  input = list(
    script = "./Notebooks/xeno/02_run_SiPSiC.Rmd",
    xeno =  pipeline$xeno_preprocess$output$xeno,
    hif_targets= hif_targets,
    genesets = msigDB_hallmarks
  ),
  output = list(
    report ="./Reports/xeno/02_run_SiPSiC/02_run_SiPSiC.html",
    sipsic_matrix ="./Reports/xeno/02_run_SiPSiC/xeno_pathwayScoresMatrix.RDS"
  )
)


pipeline[["xeno_sipsic_analysis"]] = list(
  input = list(
    script = "./Notebooks/xeno/03_SiPSiC.Rmd",
    xeno =  pipeline$xeno_preprocess$output$xeno,
    sipsic_matrix = pipeline$xenografts_run_sipsic$output$sipsic_matrix
  ),
  output = list(
    report ="./Reports/xeno/03_SiPSiC/03_SiPSiC.html",
    logFC_df = "./Reports/xeno/03_SiPSiC/xeno_sipsic_hallmarks_logFC.tsv",
    fdr_df = "./Reports/xeno/03_SiPSiC/xeno_sipsic_hallmarks_fdr_df.tsv"
  )
)

# pipeline[["xeno_sipsic_known_pathways"]] = list(
#   input = list(
#     script = "./Notebooks/xeno/06_known_pathways.Rmd",
#     xeno = "./Reports/xeno/01_preprocess/xeno.qs",
#     hallmarks_logFC_df = pipeline$xeno_sipsic_analysis$output$logFC_df,
#     hallmarks_fdr_df = pipeline$xeno_sipsic_analysis$output$fdr_df,
#     kurppa = "./input_data/pathways_from_papers/Kurppa_PMC7146079_table_s1_YAP signature.xlsx",
#     hadrek = "./input_data/pathways_from_papers/Hadrek_PMC11068778_sup_data_3.xlsx",
#     SENESCENCE = "./input_data/pathways_from_papers/FRIDMAN_SENESCENCE_UP.v2024.1.Hs.gmt",
#     maynard = "./input_data/pathways_from_papers/Maynard_PMC7484178_Table_S2.csv",
#     ferroptosis  = "./input_data/pathways_from_papers/GOBP_FERROPTOSIS.v2024.1.Hs.gmt",
#     esrra = "./input_data/pathways_from_papers/ESRRA _gene_signature.xlsx",
#     hif_targets = "./input_data/HIF_targets_Lombardi_PMC9869179.txt"
#   ),
#   output = list(
#     report ="./Reports/xeno/06_known_pathways/06_known_pathways.html"
#   )
# )

####################################### Xenogratfs cNMF ####################################################

pipeline[["xeno_cnmf_preprocess"]] = list(
  input = list(
    script = "./Notebooks/xeno/04_cnmf/01_create_data_for_cnmf.Rmd",
    xeno =  pipeline$xeno_preprocess$output$xeno
    ),
  output = list(
    report ="./Reports/xeno/04_cnmf/01_create_data_for_cnmf/01_create_data_for_cnmf.html",
    xeno_counts_filtered = "./Reports/xeno/04_cnmf/01_create_data_for_cnmf/xeno_counts_filtered.h5ad"
  ),
  params = list(
    cnmf_env = cnmf_conda_env_path
  )
)

#Note: recommended to run in background ("sbatch")
pipeline[["xeno_cnmf_run"]] = list(
  input = list(
    script = "./Notebooks/xeno/04_cnmf/02_run_cnmf/sbatch_cnmf_script.sh"
  ),
  output = list(
    cnmf_object = "Reports/xeno/04_cnmf/02_run_cnmf_1.7/models_2Kvargenes_corrected_noTPM_cnmf_obj.pckl"
  ),
  shell = substitute(paste(pipeline[[i]]$input$script))
)

pipeline[["calculate_programs"]] = list(
  input = list(
    script = "Notebooks/xeno/04_cnmf/04_calculate_programs.Rmd",
    xeno =  pipeline$xeno_preprocess$output$xeno,
    cnmf_object = "Reports/xeno/04_cnmf/02_run_cnmf_1.7/models_2Kvargenes_corrected_noTPM_cnmf_obj.pckl",
    genesets = msigDB_hallmarks
  ),
  output = list(
    report = "./Reports/xeno/04_cnmf/04_calculate_programs/04_calculate_programs.html",
    xeno_cell_usage = "./Reports/xeno/04_cnmf/04_calculate_programs/cell_usage_by_TPM.RDS",
    gep_scores = "./Reports/xeno/04_cnmf/04_calculate_programs/gep_scores.RDS"
  ),
  params = list(
    cnmf_env = cnmf_conda_env_path
  )
)


pipeline[["models_programs_analysis"]] = list(
  input = list(
    script = "./Notebooks/xeno/04_cnmf/05_models_cnmf_analysis_k5.Rmd",
    xeno =  pipeline$xeno_preprocess$output$xeno,
    xeno_cell_usage = pipeline$calculate_programs$output$xeno_cell_usage,
    gep_scores = pipeline$calculate_programs$output$gep_scores,
    hif_targets= hif_targets,
    genesets = msigDB_hallmarks
  ),
  output = list(
    report ="./Reports/xeno/04_cnmf/05_models_cnmf_analysis_k5/05_models_cnmf_analysis_k5.html"  
    )
)


####################################### patients and Bivona cNMF ####################################################
pipeline[["bivona_preprocess"]] = list(
  input = list(
    script = "./Notebooks/Bivona/01_preprocess.Rmd",
    data = "./input_data/Bivona_scRNAseq/NI04_tumor_seurat_object.RData"
  ),
  output = list(
    report ="./Reports/Bivona/01_preprocess/01_preprocess.html",
    bivona = "./Reports/Bivona/01_preprocess/bivona.qs")
)


pipeline[["patients_bivona_programs_calc"]] = list(
  input = list(
    script = "./Notebooks/patients/04_patients_cnmf_k5_from_xeno_calc.Rmd",
    lung =  pipeline$patients_preprocess$output$patients,
    bivona = pipeline$bivona_preprocess$output$bivona,
    gep_scores = pipeline$calculate_programs$output$gep_scores,
    cnmf_object = "Reports/xeno/04_cnmf/02_run_cnmf_1.7/models_2Kvargenes_corrected_noTPM_cnmf_obj.pckl"
  ),
  output = list(
    report ="./Reports/patients/04_patients_cnmf_k5_from_xeno_calc/04_patients_cnmf_k5_from_xeno_calc.html" ,
    patients_cell_usage = "./Reports/patients/04_patients_cnmf_k5_from_xeno_calc/patients_cell_usage.RDS",
    bivonas_cell_usage = "./Reports/patients/04_patients_cnmf_k5_from_xeno_calc/bivona_cell_usage.RDS"
  ),
  params = list(
    cnmf_env = cnmf_conda_env_path
  )
)

pipeline[["patients_programs_analysis"]] = list(
  input = list(
    script = "./Notebooks/patients/05_patients_cnmf_k5_from_xeno.Rmd",
    lung =  pipeline$patients_preprocess$output$patients,
    patients_cell_usage = pipeline$patients_bivona_programs_calc$output$patients_cell_usage,
    genesets = msigDB_hallmarks,
    hif_targets= hif_targets
  ),
  output = list(
    report ="./Reports/patients/05_patients_cnmf_k5_from_xeno/05_patients_cnmf_k5_from_xeno.html"  
  )
)

pipeline[["bivona_programs_analysis"]] = list(
  input = list(
    script = "./Notebooks/Bivona/02_bivona_cnmf_analysis.Rmd",
    bivona = pipeline$bivona_preprocess$output$bivona,
    bivonas_cell_usage = pipeline$patients_bivona_programs_calc$output$bivonas_cell_usage,
    genesets = msigDB_hallmarks,
    hif_targets= hif_targets
  ),
  output = list(
    report ="./Reports/Bivona/02_bivona_cnmf_analysis/02_bivona_cnmf_analysis.html"  
  )
)


####################################### Bulk cell lines - OSI ####################################################
pipeline[["bulk_cell_lines_OSI"]] = list(
  input = list(
    script = "./Notebooks/Bulk/01_cell_lines_OSI.Rmd",
    rna_counts =  "./input_data/cell_lines/OSI_bulk/OSI_bulk_cell_lines_noMTGLKI_tpm.txt",
    hif_targets= hif_targets
    ),
  output = list(
    report ="./Reports/Bulk/01_cell_lines_OSI/01_cell_lines_OSI.html"
  )
)


####################################### Bulk cell lines - OSI+ROXA June/October 23 ####################################################
# pipeline[["bulk_june23_samples_distances"]] = list(
#   input = list(
#     script = "./Notebooks/Bulk/02_bulk_osi_roxa_OCT23/01_samples_distances.Rmd",
#     rna_counts =  "./input_data/cell_lines/cell_lines_noMTGLKI_tpm.txt",
#     hif_targets= hif_targets
#   ),
#   output = list(
#     report ="./Reports/Bulk/01_cell_lines_OSI/01_cell_lines_OSI.html"
#   )
# )


####################################### Bulk cell lines - OSI+ROXA March 25 ####################################################


pipeline[["Bulk_preprocess"]] = list(
  input = list(
    script = "./Notebooks/Bulk/03_bulk_cell_lines_march25/01_preprocess.Rmd",
    rna_counts =  "./input_data/cell_lines/OSI_roxa_march25/Mar25_noMTGLKI_tpm.txt",
    sample_description = "./input_data/cell_lines/OSI_roxa_march25/RK_BIFSAMPLE.xlsx"
  ),
  output = list(
    report = "./Reports/Bulk/03_bulk_cell_lines_march25/01_preprocess/01_preprocess.html",
    bulk_Mar25_TPM_list = "./Reports/Bulk/03_bulk_cell_lines_march25/01_preprocess/bulk_Mar25_TPM_list.RDS",
    metadata_list = "./Reports/Bulk/03_bulk_cell_lines_march25/01_preprocess/metadata_list.RDS"
  )
)


pipeline[["Bulk_PCA"]] = list(
  input = list(
    script = "./Notebooks/Bulk/03_bulk_cell_lines_march25/02_PCA.Rmd",
    bulk_Mar25_TPM_list = pipeline$Bulk_preprocess$output$bulk_Mar25_TPM_list,
    metadata_list = pipeline$Bulk_preprocess$output$metadata_list
  ),
  output = list(
    report = "./Reports/Bulk/03_bulk_cell_lines_march25/02_PCA/02_PCA.html"
  )
)

pipeline[["Bulk_TPM_heatmap"]] = list(
  input = list(
    script = "./Notebooks/Bulk/03_bulk_cell_lines_march25/03_TPM_heatmap.Rmd",
    bulk_Mar25_TPM_list = pipeline$Bulk_preprocess$output$bulk_Mar25_TPM_list,
    metadata_list = pipeline$Bulk_preprocess$output$metadata_list,
    msigdb = msigdb,
    msigDB_hallmarks = msigDB_hallmarks,
    hif_targets = hif_targets,
    SASP_genes = "./input_data/Pathways/SASP_genes.txt"
  ),
  output = list(
    report = "./Reports/Bulk/03_bulk_cell_lines_march25/03_TPM_heatmap/03_TPM_heatmap.html"
  )
)

pipeline[["Bulk_GSVA"]] = list(
  input = list(
    script = "./Notebooks/Bulk/03_bulk_cell_lines_march25/04_GSVA.Rmd",
    bulk_Mar25_TPM_list = pipeline$Bulk_preprocess$output$bulk_Mar25_TPM_list,
    metadata_list = pipeline$Bulk_preprocess$output$metadata_list,
    msigdb = msigdb,
    msigDB_hallmarks = msigDB_hallmarks,
    hif_targets = hif_targets,
    SASP_genes = "./input_data/Pathways/SASP_genes.txt"
  ),
  output = list(
    report = "./Reports/Bulk/03_bulk_cell_lines_march25/04_GSVA/04_GSVA.html"
  )
)


# 
# pipeline[["bulk_mar25_HCC"]] = list(
#   input = list(
#     script = "./Notebooks/Bulk/03_bulk_cell_lines_march25/HCC_mar25_analysis.Rmd",
#     rna_counts =  "./input_data/cell_lines/OSI_roxa_march25/Mar25_noMTGLKI_tpm.txt",
#     hif_targets= hif_targets,
#     sample_description = "./input_data/cell_lines/OSI_roxa_march25/RK_BIFSAMPLE.xlsx",
#     genesets =msigDB_hallmarks,
#     SENESCENCE = "./input_data/pathways_from_papers/FRIDMAN_SENESCENCE_UP.v2024.1.Hs.gmt"
#   ),
#   output = list(
#     report ="./Reports/Bulk/03_bulk_cell_lines_march25/mar25_analysis_HCC/mar25_analysis_HCC.html",
#     up_in_persistors= "./Reports/Bulk/03_bulk_cell_lines_march25/mar25_analysis_HCC/HCC_comboVSosi_not_roxaVSctrl.txt"
#   )
# )
# 
# 
# pipeline[["bulk_mar25_H1975"]] = list(
#   input = list(
#     script = "./Notebooks/Bulk/03_bulk_cell_lines_march25/H1975_mar25_analysis.Rmd",
#     rna_counts =  "./input_data/cell_lines/OSI_roxa_march25/Mar25_noMTGLKI_tpm.txt",
#     hif_targets= hif_targets,
#     sample_description = "./input_data/cell_lines/OSI_roxa_march25/RK_BIFSAMPLE.xlsx",
#     genesets =msigDB_hallmarks,
#     SENESCENCE = "./input_data/pathways_from_papers/FRIDMAN_SENESCENCE_UP.v2024.1.Hs.gmt"
#   ),
#   output = list(
#     report ="./Reports/Bulk/03_bulk_cell_lines_march25/mar25_analysis_H1975/mar25_analysis_H1975.html",
#     up_in_persistors= "./Reports/Bulk/03_bulk_cell_lines_march25/mar25_analysis_H1975/H1975_comboVSosi_not_roxaVSctrl.txt"
#     
#   )
# )
# 
# pipeline[["march25_all_analysis"]] = list(
#   input = list(
#     script = "./Notebooks/Bulk/03_bulk_cell_lines_march25/all_cell_lines_analysis.Rmd",
#     rna_counts =  "./input_data/cell_lines/OSI_roxa_march25/Mar25_noMTGLKI_tpm.txt",
#     hif_targets= hif_targets,
#     sample_description = "./input_data/cell_lines/OSI_roxa_march25/RK_BIFSAMPLE.xlsx",
#     genesets =msigDB_hallmarks,
#     up_in_persistors_HCC = pipeline$bulk_mar25_HCC$output$up_in_persistors,
#     up_in_persistors_H1975 = pipeline$bulk_mar25_H1975$output$up_in_persistors,
#     SENESCENCE = "./input_data/pathways_from_papers/FRIDMAN_SENESCENCE_UP.v2024.1.Hs.gmt"
#     
#   ),
#   output = list(
#     report ="./Reports/Bulk/03_bulk_cell_lines_march25/all_cell_lines_analysis/all_cell_lines_analysis.html"
#   )
# )
# 



######################################## functions ###############################################

# contains_subdirectories <- function(path) {
#   items <- list.files(path, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
#   # Check which of the listed items are directories
#   # file.info() provides file/directory information, and $isdir indicates if it's a directory.
#   is_directory <- file.info(items)$isdir
#   if (any(is_directory)){
#     message("contains subdirectories")
#     user_input <- readline(prompt = "continue? y/n \n")
#     if (user_input == "y"){}
#     else {stop("aborting")}
#   }
# }



######################################## make ###############################################
library(MakefileR)
mkfile = makefile() 
all_rules = c()
mkfile = mkfile + make_rule(".Phony", "all")
for (i in 1:length(pipeline)) {
  all_rules = c(all_rules,  names(pipeline)[[i]])
}
mkfile = mkfile + make_rule("all", all_rules)

for (i in 1:length(pipeline)) {
  mkfile = mkfile +   make_comment(c("============", names(pipeline)[[i]], "============")) + # comment
    make_rule(names(pipeline)[[i]], unlist(pipeline[[i]]$output), paste("@echo  $@ is up to date") ) #define rule output
  if(is.null(pipeline[[i]]$shell)){
    mkfile = mkfile + make_rule(targets = unlist(pipeline[[i]]$output), deps = unlist(pipeline[[i]]$input), #run rscript, add & for grouped tagets
              script =paste("Rscript render.R",
                            names(pipeline)[[i]]
                            ))
  }else{
    mkfile = mkfile + make_rule(targets = unlist(pipeline[[i]]$output), deps = unlist(pipeline[[i]]$input),
                                script = eval(pipeline[[i]]$shell))
  }

}


write_makefile(makefile = mkfile,file_name = "Makefile")

