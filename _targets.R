library(targets)
library(tarchetypes)
library(future)
library(future.batchtools)

# basic check for existence of SLURM job submission command to determine if 
# running on a SLURM-enabled server
if(system('command -v sbatch') == 0) {
  plan(batchtools_slurm, template = file.path("hpc", "slurm_batchtools.tmpl"))
} else {
  plan(multisession)
}

# set packages to load
tar_option_set(
  packages = c('dplyr', 'lubridate', 'ggplot2', 'ggthemes', 'stringr', 'tidyr',
               'xtable', 'clipr', 'viridis', 'forcats'),
  deployment = 'main'
)

## load R files and workflows
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

# assemble workflow
c(
  data_targets,
  depth_bin_imputation_targets,
  dive_segmentation_targets,
  window_targets,
  pre_post_targets,
  summary_targets,
  eda_targets,
  tar_render(report_step, "results.Rmd"),
  plot_targets,
  table_targets,
  sim_plot_targets
)
