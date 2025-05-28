knitr::opts_knit$set(progress = TRUE, verbose = TRUE)
library(magrittr)
library(stringr)

pipeline = list()

####################################### Patients ####################################################


pipeline[["patients_deg"]] = list(
  input = list(
    script = "./Notebooks/patients/01_DEG.Rmd",
    lung =  "input_data/lung_cancercells_withTP_onlyPatients.rds",
    hif_targets= "./input_data/HIF_targets_Lombardi_PMC9869179.txt",
    genesets ="./input_data/h.all.v2023.2.Hs.symbols.gmt"
  ),
  output = list(
    report ="./Reports/patients/01_DEG/01_DEG.html"
  )
)

pipeline[["patients_run_sipsic"]] = list(
  input = list(
    script ="./Notebooks/patients/02_run_sipsic.Rmd",
    lung = "./input_data/lung_cancercells_withTP_onlyPatients.rds",
    hif_targets="./input_data/HIF_targets_Lombardi_PMC9869179.txt",
    genesets= "./input_data/h.all.v2023.2.Hs.symbols.gmt"
  ),
  output = list(
    report ="./Reports/patients/02_run_sipsic/02_run_sipsic.html",
    sipsic_matrix ="./Reports/patients/02_run_sipsic/patients_pathwayScoresMatrix.RDS"
  )
)

pipeline[["patients_sipsic_analysis"]] = list(
  input = list(
    script ="./Notebooks/patients/03_sipsic_analysis.Rmd",
    lung = "./input_data/lung_cancercells_withTP_onlyPatients.rds",
    sipsic_matrix = pipeline$patients_run_sipsic$output$sipsic_matrix
  ),
  output = list(
    report ="./Reports/patients/03_sipsic_analysis/03_sipsic_analysis.html"
  )
)

####################################### Xenogratfs ####################################################
pipeline[["xenografts_deg"]] = list(
  input = list(
    script = "./Notebooks/xeno/01_DEG.Rmd",
    xeno =  "input_data/xeno.qs",
    hif_targets= "./input_data/HIF_targets_Lombardi_PMC9869179.txt",
    genesets ="./input_data/h.all.v2023.2.Hs.symbols.gmt"
  ),
  output = list(
    report ="./Reports/xeno/01_DEG/01_DEG.html"
  )
)

pipeline[["xenografts_run_sipsic"]] = list(
  input = list(
    script = "./Notebooks/xeno/02_run_SiPSiC.Rmd",
    xeno =  "input_data/xeno.qs",
    hif_targets= "./input_data/HIF_targets_Lombardi_PMC9869179.txt",
    genesets ="./input_data/h.all.v2023.2.Hs.symbols.gmt"
  ),
  output = list(
    report ="./Reports/xeno/02_run_SiPSiC/02_run_SiPSiC.html",
    sipsic_matrix ="./Reports/xeno/02_run_SiPSiC/xeno_pathwayScoresMatrix.RDS"
  )
)


pipeline[["xeno_sipsic_analysis"]] = list(
  input = list(
    script = "./Notebooks/xeno/03_SiPSiC.Rmd",
    xeno =  "input_data/xeno.qs",
    sipsic_matrix = pipeline$xenografts_run_sipsic$output$sipsic_matrix
  ),
  output = list(
    report ="./Reports/xeno/03_SiPSiC/03_SiPSiC.html"
  )
)
####################################### Bulk cell lines ####################################################


pipeline[["bulk_mar25_HCC"]] = list(
  input = list(
    script = "./Notebooks/Bulk/bulk_cell_lines_march25/HCC_mar25_analysis.Rmd",
    rna_counts =  "./input_data/osiRoxa_bulk/Mar25/gene_count.xls",
    hif_targets= "./input_data/HIF_targets_Lombardi_PMC9869179.txt",
    sample_description = "./input_data/osiRoxa_bulk/Mar25/RK_BIFSAMPLE.xlsx",
    genesets ="./input_data/h.all.v2023.2.Hs.symbols.gmt"
  ),
  output = list(
    report ="./Reports/Bulk/bulk_cell_lines_march25/mar25_analysis_HCC.html"
  )
)


pipeline[["bulk_mar25_H1975"]] = list(
  input = list(
    script = "./Notebooks/Bulk/bulk_cell_lines_march25/H1975_mar25_analysis.Rmd",
    rna_counts =  "./input_data/osiRoxa_bulk/Mar25/gene_count.xls",
    hif_targets= "./input_data/HIF_targets_Lombardi_PMC9869179.txt",
    sample_description = "./input_data/osiRoxa_bulk/Mar25/RK_BIFSAMPLE.xlsx",
    genesets ="./input_data/h.all.v2023.2.Hs.symbols.gmt"
  ),
  output = list(
    report ="./Reports/Bulk/bulk_cell_lines_march25/mar25_analysis_H1975.html"
  )
)
######################################## functions ###############################################
#get input/output from pipeline with script name
get_input <- function(script) {
  for (i in seq_along(pipeline)) {
    item <- pipeline[[i]]
    if (is.list(item) && "input" %in% names(item) && "script" %in% names(item$input)) {
      if (item$input$script == script) {
        return(pipeline[[i]]$input)
      }
    }
  }
}

get_output <- function(script) {
  for (i in seq_along(pipeline)) {
    item <- pipeline[[i]]
    if (is.list(item) && "input" %in% names(item) && "script" %in% names(item$input)) {
      if (item$input$script == script) {
        return(pipeline[[i]]$output)
      }
    }
  }
}

get_label <- function(script) {
  for (i in seq_along(pipeline)) {
    item <- pipeline[[i]]
    if (is.list(item) && "input" %in% names(item) && "script" %in% names(item$input)) {
      if (item$input$script == script) {
        return(names(pipeline)[[i]])
      }
    }
  }
}

get_current_path <- function() {
  # Get the full path of the current file
  full_path <- rstudioapi::getSourceEditorContext()$path
  
  # Get the full path of the project directory
  project_path <- rstudioapi::getActiveProject()
  
  # Normalize the paths to ensure consistency
  normalized_full_path <- normalizePath(full_path,winslash = "/")
  normalized_project_path <- normalizePath(project_path,winslash = "/")
  
  # Remove the project directory path from the full file path
  # and replace it with "./" to create the relative path
  relative_path <- gsub(paste0("^", normalized_project_path, "/"), "./", normalized_full_path)
  
  # Print the relative path
  return(relative_path)
}

my_render <- function(notebook_path , set_params = list()){
    knitr::opts_knit$set(progress = TRUE, verbose = TRUE)
    report = get_output(script = notebook_path)$report
    input = get_input(notebook_path)
    output = get_output(notebook_path)
    if (dir.exists(dirname(report))) {
      unlink(dirname(report),recursive = T)
    }
    set_params[["data_out_dir"]] = dirname(report) %s+% "/"
    message("Rendering to:")
    message(dirname(report))
    
    rmarkdown::render(
      input = notebook_path,
      output_format = "html_document",
      output_file = report,
      knit_root_dir = getwd(),
      output_dir = dirname(report),
      params = set_params
    )


}

#################################### makepipe ####################################################################


library(makepipe)
makepipe::reset_pipeline()
makepipe_pipe <- get_pipeline()


script = pipeline[[1]]$input$script
make_with_recipe(
  recipe = my_render(notebook_path =pipeline[[1]]$input$script),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = get_label(script),build = F
)



script = pipeline[[2]]$input$script
make_with_recipe(
  recipe = my_render(notebook_path =pipeline[[2]]$input$script),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = get_label(script),build = F
)


script = pipeline[[3]]$input$script
make_with_recipe(
  recipe = my_render(notebook_path =pipeline[[3]]$input$script),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = get_label(script),build = F
)

script = pipeline[[4]]$input$script
make_with_recipe(
  recipe = my_render(notebook_path =pipeline[[4]]$input$script),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = get_label(script),build = F
)

script = pipeline[[5]]$input$script
make_with_recipe(
  recipe = my_render(notebook_path =pipeline[[5]]$input$script),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = get_label(script),build = F
)

script = pipeline[[6]]$input$script
make_with_recipe(
  recipe = my_render(notebook_path =pipeline[[6]]$input$script),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = get_label(script),build = F
)