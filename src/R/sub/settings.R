#' Settings
#' create : 2020/10/15


dir.import_data <- "../../mail/200817/data/res"
fn.import_data  <- "col_info_added_received_data.xlsx"

dir.RData <- "../data"
fn.imported_data <- "imported_data.RData"

library(grDevices)

# Setting for histogram -----------------
#' functions for argument "binwidth=..."
bw = function(x, f) ceiling((max(x) - min(x)) / f(x))
FD = function(x) bw(x, nclass.FD)
scott = function(x) bw(x, nclass.scott)
sturges = function(x) bw(x, nclass.Sturges)
# endrant --------------------------------