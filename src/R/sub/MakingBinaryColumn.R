#' (one phrase discription what this function does.)
#' @import magritrr
#' @import tidyr
#' @import plyr
#' @import dplyr
#' @import tibble
#' 
#' @parameter df A data.frame object.
#' @parameter columnnumber A numeric vector indexing target column.
#' @parameter cutoff cutoff value you choose for variable (e.g. 71, you can use "1" for 0 or 1 value)
#' 
#' @export

MakingBinaryColumn <- function(
  df,columnnumber, cutoff
  ) {
  NewColumnName<-sprintf(as.character(columnnumber),cutoff)
  for (i in 1:nrow(df)) {
    if (is.na(df[,columnnumber][i])==TRUE) {df$NewColumnName[i] <-"NA"} else if (df[,columnnumber][i] < cutoff) {df$NewColumnName[i] <-1
    } else {
      df$NewColumnName[i] <-0
    }
  }
  df$NewColumnName
}