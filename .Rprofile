set.seed(3320)
library(Matrix)
library(stringi)
library(rlang)
library(Seurat)
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
# conflict_prefer(name = "list", winner = "base")
library(rmarkdown)
library(magrittr)
# source("./build_funs.R")


# my_render <- function(script, report,input,output, input_names, output_names ) {
#   input = strsplit(x = input,split = " ") %>% as.list()
#   names(input) = input_names
#   output = strsplit(x = output,split = " ") %>% as.list()
#   names(output) = output_names
#   rmarkdown::render(script , output_dir = dirname(report),knit_root_dir = getwd(), 
#                     params = list(data_out_dir =dirname(report)))
#   # Notebooks/xeno/01_DEG.Rmd input_data/xeno.qs input_data/h.all.v2023.2.Hs.symbols.gmt input_data/HIF_targets_Lombardi_PMC9869179.txt
# }
