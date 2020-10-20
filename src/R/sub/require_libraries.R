# options(BioC_mirror="https://bioconductor.org/")
# source("https://bioconductor.org/biocLite.R")


packages.in.CRAN <- c(
  "tidyverse", "ggplot2", "reshape2","readxl", "dbplyr", "tableone",
  "pROC","Matching","survey"
  )

for(i in 1:length(packages.in.CRAN)){
  if (!requireNamespace(packages.in.CRAN[i], quietly = TRUE)) install.packages(packages.in.CRAN[i])
  eval(
    parse(text=sprintf("require(%s)", packages.in.CRAN[i]))
  )
}


# if(!require(ExploratoryDataAnalysis)){
#   devtools::install_github("mrmtshmp/ExploratoryDataAnalysis")
# }


if(Bibtex){
  write(toBibtex(citation()),file="CRAN")
  for(i in 1:length(packages.in.CRAN)){
    write(toBibtex(citation(packages.in.CRAN[i])),file=sprintf("./src/biblio/%s%s.bib",packages.in.CRAN[i],"_CRAN"))
  }
}

