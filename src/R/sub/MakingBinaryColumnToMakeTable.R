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
#' @parameter selected_var column numbers you choose to compare between two groups (e.g. c(2,4,6,8,9,16))
#' 
#' @export

MakingBinaryColumnToMakeTable <- function(
  df,
  columnnumber,
  cutoff,
  selected_var
  ) {
  NewColumnName<-sprintf(as.character(columnnumber),cutoff)
  for (i in 1:nrow(df)) {
    if (is.na(df[,columnnumber][i])==TRUE) {
      df$NewColumnName[i] <-"NA"
    } else if (df[,columnnumber][i] < cutoff) {
      df$NewColumnName[i] <-1
    } else {
      df$NewColumnName[i] <-0
    }
  }
  NEW_pos<-
    sprintf(as.character(columnnumber),cutoff,"_pos") 
  NEW_neg<-
    sprintf(as.character(columnnumber),cutoff,"_neg") 
  NEW_pos<-
    df %>% subset(NewColumnName==1)
  NEW_neg<-
    df %>% subset(NewColumnName==0)
  NEW_result<-
    c("Name","positive_q1","positive_median","positive_q3","negative_q1","negative_median","negative_q3","p.value", "n.pos", "n.neg") 
  
  for (i in selected_var){
    i<-c((colnames(df)[i]),getElement(quantile(NEW_pos[,i],1/4,na.rm=TRUE),"25%"),getElement(quantile(NEW_pos[,i],2/4,na.rm=TRUE),"50%"),getElement(quantile(NEW_pos[,i],3/4,na.rm=TRUE),"75%"),getElement(quantile(NEW_neg[,i],1/4,na.rm=TRUE),"25%"),getElement(quantile(NEW_neg[,i],2/4,na.rm=TRUE),"50%"),getElement(quantile(NEW_neg[,i],3/4,na.rm=TRUE),"75%"),(getElement(wilcox.test(NEW_pos[,i],NEW_neg[,i],exact=FALSE),"p.value")),length(na.omit(NEW_pos[,i])),length(na.omit(NEW_neg[,i])))
    NEW_result<-rbind(NEW_result,i)
  }
  
  NEW_result<-data.frame(NEW_result)
  NEW_result
}