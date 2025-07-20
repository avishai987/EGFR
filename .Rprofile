set.seed(3320)
check_for_uncommitted_changes <- function(repo_path = ".") {
  
  # Check if a Git repository exists at the specified path
  if (!dir.exists(file.path(repo_path, ".git"))) {
    message("⚠️ No Git repository found at the specified path: ", repo_path)
    return(invisible(FALSE))
  }
  
  # Run 'git status --porcelain' to check for uncommitted changes
  # '--porcelain' provides a stable, easy-to-parse output.
  # If there are changes, the output will have lines; otherwise, it's empty.
  git_status_output <- system("git status --porcelain", intern = TRUE, ignore.stderr = TRUE)
  
  # Check if there's any output (meaning uncommitted changes exist)
  if (length(git_status_output) > 0) {
    message("🚨 You have uncommitted changes in '", basename(normalizePath(repo_path)), "'!")
    message("  Please run 'git add' and 'git commit' to save your work.")
    message("\n  Details of changes:")
    cat(paste0("  ", git_status_output, collapse = "\n"))
    message("\n")
    return(TRUE)
  } else {
    message("✅ No uncommitted changes. Great job!")
    return(FALSE)
  }
}


check_for_uncommitted_changes()
my_library_folder = "/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/libs"
.libPaths(c(my_library_folder, .libPaths())) # make my library folder as default folder (default folder is not writable)
ulimit::memory_limit(30000) # Limit ram to avoid crash
library(igraph, lib.loc = "/usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/r-igraph-1.4.2-lxqmjthw2lo45ggzwhg3x47ehnezig22/rlib/R/library")
library(Matrix)
library(stringi)
library(rlang)
library(Seurat,lib.loc = "../libs/seurat_4.0.1/")
library(ggplot2)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(devtools)
library(SourceFromGithub)
library("plyr")
library(dplyr)
library(stringi)
library(conflicted)
library(ggplotify)
library(print.rmd.tabs)
conflict_prefer(name = "intersect", winner = "base")
conflict_prefer(name = "select", winner = "dplyr")
conflict_prefer(name = "filter", winner = "dplyr")
conflict_prefer(name = "mutate", winner = "dplyr")
conflict_prefer(name = "arrange", winner = "dplyr")
conflicts_prefer(dplyr::summarize)

conflicts_prefer(dplyr::desc)
library(rmarkdown)
library(magrittr)



