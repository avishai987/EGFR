options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("MakefileR")
devtools::install_github("avishai987/print.rmd.tabs")
devtools::install_github("avishai987/SourceFromGithub")
devtools::install_github("montilab/hypeR@e407bf1", dependencies = F) # cannot install with conda due to https://github.com/montilab/hypeR/issues/58
devtools::install_github("iholzleitner/facefuns@5a7df7d", upgrade = F)
