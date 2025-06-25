library(reticulate)
use_condaenv("bioenv", conda = "/opt/conda/condabin/conda", required = TRUE)
# py_list_packages()
reticulate::repl_python()
