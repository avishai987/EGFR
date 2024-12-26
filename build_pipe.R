knitr::opts_knit$set(progress = TRUE, verbose = TRUE)
library(magrittr)
library(stringr)
get_report_target_path <- function(notebook_path) {
  filename = gsub(x = notebook_path,pattern = "\\.Rmd",replacement = ".html")
  dir_name = str_remove(string = notebook_path, pattern = "./Notebooks/") %>% str_remove(pattern = ".Rmd")
  new_dir_path = paste0("./Reports/", dir_name, "/")
  report_path =  paste0(new_dir_path,basename(filename))
  return(report_path)
}

get_data_target_path <- function(notebook_path,data_name) {
  filename = gsub(x = notebook_path,pattern = "\\.Rmd",replacement = "")
  dir_name = str_remove(string = filename, pattern = "./Notebooks/") 
  new_dir_path = paste0("./Reports/", dir_name, "/")
  data_path =  paste0(new_dir_path,data_name)
  return(data_path)
}


get_targets_and_depens <- function(notebook_path) {
  yml_metadata <- rmarkdown::yaml_front_matter(notebook_path)
  dependencies = yml_metadata$params$input_data
  if (!is.null(dependencies)) {
    dependencies = dependencies  %>% lazyeval::lazy_eval() %>% base::unname()
  }
  
  targets = yml_metadata$params$output_data
  if (!is.null(targets)) {
    targets = targets %>% lazyeval::lazy_eval() %>%
      get_data_target_path(notebook_path = notebook_path)
  }
  
  report_path = get_report_target_path(notebook_path)
  targets = c(targets,report_path)
  return(list(targets = targets, dependencies = dependencies))
}

my_render <- function(notebook_path , set_params = list()) 
{
  filename = gsub(x = notebook_path,pattern = "\\.Rmd",replacement = ".html")
  dir_name = str_remove(string = notebook_path, pattern =  "./Notebooks/") %>% str_remove(pattern = ".Rmd")
  new_dir_path = paste0("./Reports/", dir_name, "/")
  if(dir.exists(new_dir_path)) {
    unlink(new_dir_path, recursive = T)
  }
  dir.create(path = new_dir_path, recursive = T)
  
  set_params[["data_out_dir"]] = new_dir_path
  message("Rendering to:")
  message(new_dir_path)
  
  rmarkdown::render(
    input = notebook_path,
    output_format = "html_document",
    output_file = filename,
    knit_root_dir = getwd(),
    output_dir = new_dir_path,
    params = set_params)
  
  
}

add_segment <- function(recipe, targets, dependencies = NULL, packages = NULL, 
                        envir = new.env(parent = parent.frame()), quiet = getOption("makepipe.quiet"), 
                        force = FALSE, label = NULL, note = NULL) {
  library(makepipe)
  recipe <- substitute(recipe)
  pipeline <- get_pipeline()
  if (is.null(pipeline)) {
    pipeline <- Pipeline$new()
    set_pipeline(pipeline)
  }
  segment <- pipeline$add_recipe_segment(recipe, targets, 
                                         dependencies, packages, envir, force)
  makepipe:::add_note_and_label(pipeline, segment, label, note)
}

create_pipe <- function() {
  library(makepipe)
  makepipe::reset_pipeline()
  pipe <- get_pipeline()
  # add a segment to the pipeline
  add_segment(
    recipe = my_render(notebook_path = "./Notebooks/test.Rmd",set_params = list(my_print = "ooo")),
    targets = "./Reports/test/test.html",
    dependencies = "./Notebooks/test.Rmd")
}

