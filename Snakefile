rule models_deg:
    input: 
       script = "Notebooks/xeno/01_DEG.Rmd",
       xeno = "input_data/xeno.qs",
       genesets = "input_data/h.all.v2023.2.Hs.symbols.gmt",
       hif_targets= "input_data/HIF_targets_Lombardi_PMC9869179.txt",

    output:
      report = "Reports/xeno/01_DEG/01_DEG.html"
    shell:
      """
      Rscript -e \"my_render ( '{input.script}', \"
      """

################ Patients ######################################
rule patients_deg:
    input: 
      script = "Notebooks/patients/01_DEG.Rmd",
      lung =  "input_data/lung_cancercells_withTP_onlyPatients.rds",
      hif_targets= "input_data/HIF_targets_Lombardi_PMC9869179.txt",
      genesets=  "input_data/h.all.v2023.2.Hs.symbols.gmt"
    output:
      report = "Reports/patients/01_DEG/01_DEG.html"
    shell:
      """
      Rscript -e \"env = new.env(); env$input = list(lung = '{input.lung}',hif_targets =  '{input.hif_targets}', genesets = '{input.genesets}') ; rmarkdown::render('{input.script}',envir=env, output_dir = dirname('{output.report}'),knit_root_dir = getwd(),params = list(data_out_dir =dirname('{output.report}'))) \"
      """
      
      
rule patients_run_sipsic:
    input: 
      script = "Notebooks/patients/02_run_sipsic.Rmd",
      lung =  "input_data/lung_cancercells_withTP_onlyPatients.rds",
      hif_targets= "input_data/HIF_targets_Lombardi_PMC9869179.txt",
      genesets=  "input_data/h.all.v2023.2.Hs.symbols.gmt"
    output:
      report = "Reports/patients/02_run_sipsic/02_run_sipsic.html",
      sipsic_matrix = "Reports/patients/02_run_sipsic/patients_pathwayScoresMatrix.RDS"
    shell:
      """
      Rscript -e \"env = new.env(); env$input = list(lung = '{input.lung}',hif_targets =  '{input.hif_targets}', genesets = '{input.genesets}') ;\
env$output = list(sipsic_matrix = '{output.sipsic_matrix}'); \
rmarkdown::render('{input.script}',envir=env, output_dir = dirname('{output.report}'),knit_root_dir = getwd(), \
params = list(data_out_dir =dirname('{output.report}'))) \"
      """
      
rule patients_sipsic_analysis:
    input: 
      script = "Notebooks/patients/03_sipsic_analysis.Rmd",
      lung =  "input_data/lung_cancercells_withTP_onlyPatients.rds",
      sipsic_matrix = rules.patients_run_sipsic.output.sipsic_matrix
    output:
      report = "Reports/patients/03_sipsic_analysis/03_sipsic_analysis.html",
    shell:
      """
      Rscript -e \"env = new.env(); env$input = list(lung = '{input.lung}',sipsic_matrix =  '{input.sipsic_matrix}') ;\
rmarkdown::render('{input.script}',envir=env, output_dir = dirname('{output.report}'),knit_root_dir = getwd(), \
params = list(data_out_dir =dirname('{output.report}'))) \"
      """


rule all:
  input:
    report1 = rules.patients_deg.output.report,
    report2 = rules.patients_run_sipsic.output.report,
    report3 = rules.patients_sipsic_analysis.output.report

    
    
