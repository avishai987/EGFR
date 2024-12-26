
build_pipe <- function() {
  job::job({  
    source(".Rprofile")
    source("build_pipe.R")
    create_pipe()
    pipe <- get_pipeline()
    pipe$build()
  },import = NULL,packages = NULL)
  
}

plot_pipe <- function() {
  source(".Rprofile")
  source("build_pipe.R")
  create_pipe()
  makepipe::show_pipeline(as = "visnetwork")
}