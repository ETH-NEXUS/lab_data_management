## Path where the report file should go
output_dir = "C:/Users/mprummer/OneDrive - ETH Zurich/00_NEXUS/08_groups/02_st/Projects/2025_st_internal_LDM/"

## Name of the report file
output_file = "241126_nexusAD_SW620_report" # or NULL: output_file = rmd_file
## Format of the report file 
#output_format = "html_document" # or "pdf_document" or "word_document"
output_format = "html" # or "pdf" or "word"
## Path to SLmisc.R
path_SLmisc = "C:/Users/mprummer/Documents/_Michael/nexus_github/SLmisc/R/"
## Path to the analysis scripts
path_Rscript = "C:/Users/mprummer/OneDrive - ETH Zurich/00_NEXUS/08_groups/02_st/Projects/2025_st_internal_LDM/"

params = list(
  project = "2024_snl_pruschy_radiation",
  screen = "241126_Screen_NexusApproved_SW620",
  hts_type = "selectivity", # c("single", "selectivity", "doseresponse")
  condi_yes = "irradiated",
  condi_no = "not irradiated",
  geom_cor = "median polish", # c("no", "median polish", "local polynomial fit", "loess")
  act_cut = "log10(1.5)", # must be a string
  fdr_cut = 0.01, # [0, 1]
  select_cut = 0.3, # [0, oo]
  select_cut_yes = "2 * act_cut", # target act, must be a string
  select_cut_no = "act_cut", # counter act, must be a string
  # path to library description "CHEMICAL INFO"
  path_lib = "C:/Users/mprummer/Downloads/241126_Screen_NexusApproved_SW620_Lum.csv",
  # path to measurement data "MAIN INFO"
  path_data = "C:/Users/mprummer/Downloads/241126_Screen_NexusApproved_SW620_Lum (1).csv",
  # path to experiment description "EXPERIMENT DATA"
  path_meta = "C:/Users/mprummer/OneDrive - ETH Zurich/00_NEXUS/08_groups/02_st/Projects/_old projects/2024_snl_pruschy_radiation/241126_Screen_NexusApproved_MC38_SW620/SW620/experiment_data.csv",
  # output path for all result files
  path_output = "C:/Users/mprummer/OneDrive - ETH Zurich/00_NEXUS/08_groups/02_st/Projects/_old projects/2024_snl_pruschy_radiation/241126_Screen_NexusApproved_MC38_SW620/SW620/",
  # path to SLmisc.R
  path_SLmisc = path_SLmisc
)


## install / load required packages
# if (!requireNamespace("SLmisc", quietly = TRUE)) {
#   remotes::install_github("ETH-NEXUS/SLmisc")
# }
# library(SLmisc)
source(paste0(path_SLmisc, "SLmisc.R") )
lby = c("RColorBrewer", "ggplot2", "gplots", "ggrepel", "quarto", "cowplot",
        "ggpubr", "ggExtra", "dplyr", "reshape2", "lubridate",
        "reshape2", "knitr", "pheatmap", "tidyr", "locfit", "robust","MASS")
for(ii in seq_len(length(lby))) {
  if (!requireNamespace(lby[ii], quietly = TRUE)) {
    install.packages(lby[ii])
  }
  resp = require(lby[ii], character.only=T, warn.conflicts=F, quietly=T, lib.loc=.libPaths()[1])
  if(!resp) stop("Could not load one or more required packages")
}
rm(resp, lby)



## suggestion for options:
if(params$hts_type == "single"){
  rmd_file = path_Rscript %&% "single.qmd"
}

if(params$hts_type == "selectivity"){
  rmd_file = path_Rscript %&% "selectivity.qmd"
}

if(params$hts_type == "doseresponse"){
  rmd_file = path_Rscript %&% "doseresponse.qmd"
}
## end suggestion


## render the qmd file, produce report file.
quarto::quarto_render(input = rmd_file, output_file = output_file, 
                  output_format = output_format, 
                  execute_params = params, quiet = F)
if(path_Rscript != output_dir) {
  file.copy(path_Rscript %&% output_file %&% "." %&% output_format, output_dir %&% output_file %&% "." %&% output_format)
  file.remove(path_Rscript %&% output_file %&% "." %&% output_format)
}
write.table(reshape2::melt(unlist(params)), output_dir %&% "params.tsv", sep = "\t")
