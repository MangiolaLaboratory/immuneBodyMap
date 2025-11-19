library(targets)
library(tidyverse)
library(cli)
library(glue)

# cell type:
# [1] "b naive"         "cd14 mono"       "cd16 mono"       "cd4 naive"
# [5] "cd4 tcm"         "cd4 th17 em"     "cd8 naive"       "cd8 tcm"
# [9] "cd8 tem"         "t cd4"           "t cd8"           "tgd"
# [13] "treg"            "monocytic"       "b memory"        "cd4 th1/th17 em"
# [17] "cytotoxic"       "cd4 th2 em"      "progenitor"      "mait"
# [21] "nk"              "t"               "macrophage"      "cd4 fh em"
# [25] "cd4 tem"         "cd4 th1 em"      "cdc"             "b"
# [29] "dc"              "plasma"          "granulocyte"     "erythrocyte"
# [33] "pdc"             "ilc"             "neuron"          "glial"
# [37] "pericyte"        "endothelial"     "immune"          "blood"
# [41] "muscle"          "stromal"         "mesothelial"     "epithelial"
# [45] "liver"           "mast"            "nkt"             "renal"
# [49] "endocrine"       "reproductive"    "secretory"       "fat"
# [53] "pneumocyte"      "myoepithelial"   "sensory"         "lens"
# [57] "epidermal"       "cartilage"       "bone"

# cell_type_list <- c(
#   "plasma", 
#   "nk",
#   "cd4 naive",
#   "cd8 naive",
#   "cd4 th1/th17 em",
#   "cd4 th17 em",
#   "mait",
#   "cytotoxic",
#   "tgd",
#   "cd4 fh em",
#   "cd4 th2 em",
#   "cd8 tem",
#   "cd4 tcm",
#   "b memory",
#   "b naive",
#   "macrophage",
#   "cd8 tcm",
#   "treg"
# )

cell_type_list <- 'epithelial'

# cell_type_list <- c(
#   'stromal',
#   'muscle',
#   'secretory',
#   'pericyte',
#   'endothelial'
# )

# Store results
done_list <- list()
error_list <- list()

root_path = '/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE'
script_path <- '/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/TAR_SCRIPTS/dynamic_tar_script.R'

cli_progress_bar("Processing cell types", total = length(cell_type_list))

for (cur_ct in cell_type_list) {
  
  setwd(root_path)
  
  tryCatch({
    
    cli_alert_info(glue('Working on {cur_ct}...'))
    
    cur_path = glue::glue('{root_path}/V1_{cur_ct %>% make.names}')
    if(!dir.exists(cur_path)) dir.create(cur_path)
    if(dir.exists(cur_path)) cli_alert_success(glue('{cur_ct} directory created successfully!'))
    
    setwd(cur_path)
    source(script_path)
    tar_make(hypothesis_sex_differ_at_decades)
    
    cli_alert_success(glue('{cur_ct} ALL done!'))
    done_list <- append(done_list, cur_ct)
    
  }, error = function(e) {
    cli_alert_danger(glue("Error in {cur_ct}: {e$message}"))
    error_list[[cur_ct]] <<- e$message
  })
  
  # Progress update regardless of success or error
  cli_progress_update()
  
}

cli_progress_done()