Sys.setenv(TMPDIR="./tmp")
configure.vars="TMPDIR=./tmp"
list.of.packages <- c("stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos="http://cran.r-project.org")
