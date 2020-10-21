# ---
# title: "9Aug2020"
# author: "Saeko Fukui"
# date: "8/9/2020"
# output: html_document
# ---

#' "ATT" association of *C3>71 at 2w from LDLT (varname:C3_2w_71)* with *mortality90*.
#' The propensity score was calculated as predicted logit of
#' high level C3 at 2w after LDLT (*C3_2w_71==1*) .
#' The predictors for the logit were *PCT_pre*, *infection14* and *IVIg*.
#' 
#' **Definition of the variables:**
#' *C3_2w_71* <-  ## The primary endpoint. 
#'   MakingBinaryColumn(
#'     LDLT_Comp_completed,
#'     which("C3.2w" == colnames(LDLT_Comp_completed)),
#'     cutoff = 71
#'   )
#'

fml.ps_model <- as.factor(C3_2w_71) ~
  PCT_pre + infection14 + IVIg

##```{r #original data to use}

source(
  "require_libraries.R"
  )

source(
  "settings.R"
  )

##```
  
##```{r #original data to use}

LDLT_Comp_completed <- 
  read_excel(
    "original data/LDLT_Comp_completed.xlsx",
    sheet = "Sheet1", na = "NA"
    ) %>%
  data.frame()
# LDLT_Comp_completed<-data.frame(LDLT_Comp_completed)
# colnames(LDLT_Comp_completed)
column_for_Fisher<-c(3,4,5,7,8,14,15,25,34,35,36,37,38,39,40,163,166,168,170,171,172,174,176,178,180,182,184,186,190,197,198)
##```

##```{r for Table 1 and Table 2}

# Load source files of functions.
dir.Rsources <- './prg/R'
list.fn.Rsources <- list.files(dir.Rsources)

for(i in 1:length(list.fn.Rsources))
  source(
    sprintf(
      '%s/%s',dir.Rsources,list.fn.Rsources[i]
    )
  )

selected_var<-c(2,6,9, 10,11,12,13,16:19,41:160)

write.table(MakingBinaryColumnToMakeTable(LDLT_Comp_completed, 198, 1, selected_var), "results/90_mortality_wilcoxon.csv",quote=F,col.names=F, row.names=F, append=F)
write.table(MakingBinaryColumnToMakeTable(LDLT_Comp_completed, 180, 1, selected_var), "results/7_infection_wilcoxon.csv",quote=F,col.names=F, row.names=F, append=F)
write.table(MakingBinaryColumnToMakeTable(LDLT_Comp_completed, 111, 71, selected_var), "results/C3_2w_71_wilcoxon.csv",quote=F,col.names=F, row.names=F, append=F)

LDLT_Comp_completed$C3_2w_71 <- 
  MakingBinaryColumn(
    LDLT_Comp_completed,
    111,
    71
    )
LDLT_Comp_completed <-
  LDLT_Comp_completed %>%
  dplyr::filter(!is.na(C3.2w))

#```
##```{r #compare binary results for selected variables for Table 1 and Table 2}

column_for_Fisher <-
  c(3,5,7,8,14,15,25,34,35,37,38,39,40,56,163,166,168,170,171,172,174,176,178,180,182,184,186,197,198,199)

## Covariates
vars <- colnames(LDLT_Comp_completed)[column_for_Fisher]

## Construct a table
tabUnmatched <-
  CreateTableOne(
    vars = vars, 
    strata = "C3_2w_71", 
    data = LDLT_Comp_completed, 
    test = FALSE)
## Show table with SMD
sink(
  "results/TableOne.txt"
  )
print(tabUnmatched, smd = TRUE)
sink()



# write.table(
#   FisherAll(
#     LDLT_Comp_completed, column_for_Fisher
#     ),
#   "results/Fisher_for_Table1&2.csv",
#   quote=F,col.names=F, row.names=F, append=F
#   )

#```
##```{r}
for (
  i in 1:nrow(LDLT_Comp_completed)
  ) {
  if(LDLT_Comp_completed$Days_of_lifetime[i] <= 90) {
    LDLT_Comp_completed$Time_mortality90[i] <-LDLT_Comp_completed$Days_of_lifetime[i]
      } else {
        LDLT_Comp_completed$Time_mortality90[i] <-90
      }
  }
#```
##```{r}

# Propensity score model ---------------

propensityScoreModel <- 
  glm(
    fml.ps_model,
    family  = binomial(link = "logit"),
    data    = LDLT_Comp_completed
    )

LDLT_Comp_completed.propensityScores <-
  LDLT_Comp_completed %>%
  rownames_to_column("seqid") %>%
  left_join(
    predict(
      propensityScoreModel,
      type="response",
      na.action = na.exclude
      ) %>%
      data.frame() %>%
      rownames_to_column("seqid"),
    by="seqid"
    ) %>%
  dplyr::rename(
    "propensity_score" = "."
    ) %>%
  mutate( # Add weight columns calculated from the propensity score.
    C3_2w_71 = as.numeric(C3_2w_71)
    )

LDLT_Comp_completed.propensityScores_IPW <-
  IPW_weights(
    treatment = LDLT_Comp_completed.propensityScores$C3_2w_71,
    propensity_score = LDLT_Comp_completed.propensityScores$propensity_score,
    dat = LDLT_Comp_completed.propensityScores
    )

res.roc.propensity_score <- roc(
  response = 
    as.factor(LDLT_Comp_completed.propensityScores_IPW$C3_2w_71),
  predictor = 
    LDLT_Comp_completed.propensityScores_IPW$propensity_score
  )

# Plot distribution of propensity score and weighted counts. ----------------
#'
ggdata.propensityScores <-
  ggplot(
    data =
      LDLT_Comp_completed.propensityScores_IPW,
    aes(
      x = propensity_score
      )
    )
ggdata.propensityScores.weighted_count <-
  ggplot(
    data =
      LDLT_Comp_completed.propensityScores_IPW %>%
      pivot_longer(
        cols =
          c(starts_with("w_at")),
        values_to = "weight",
        names_to = "target_pop"
      ) %>%
      dplyr::filter(target_pop=="w_att"),
    aes(
      x = propensity_score,
      weight = weight
    )
  )

pdf(
  file = "results/IPWcount.pdf",width = 21
  )
plot(
  ggdata.propensityScores + 
    geom_density(
      aes(
        fill = 
          as.factor(C3_2w_71)
        ),
      bw="SJ",
#      binwidth = FD,
      alpha=0.5,
      position="identity"
      ) + 
    theme_bw()
  )
plot(
  ggdata.propensityScores.weighted_count + 
    geom_density(
      aes(fill=as.factor(C3_2w_71)),
      bw="SJ",
#      binwidth = FD,
      alpha=0.5#,
#      position="dodge"
      ) + 
    facet_grid(~target_pop) + 
    theme_bw()
    )
plot(
  res.roc.propensity_score,
  print.thres=TRUE 
  )
legend(
  x = 0.6, y=0.5,cex = 0.7, 
  # lwd = c(2,2,0), lty = 1:2,
  legend = c(
    sprintf(
      "AUC = %s (0.95CI: %s, %s)",
      round(auc(res.roc.propensity_score),3),
      round(ci(auc(res.roc.propensity_score))[1],3),
      round(ci(auc(res.roc.propensity_score))[3],3)
    )
  ),
  bty = "n"
)
dev.off()


# Balancing assessment ----------------------------------------------------

#' The standardized mean differences between the two patients' groups of *low serum C3* and *high serum C3* 
#' in each of covariates those which 
#' The tableone package (tableone_CRAN.bib)


#' Reference:
#' https://cran.r-project.org/web/packages/tableone/vignettes/smd.html
#' (Many lines were snipped from above website [accessed:2020/08/31])

res.svydesign.w_att <- 
  survey::svydesign(
    ids = ~ 1, 
    data = LDLT_Comp_completed.propensityScores_IPW[!is.na(LDLT_Comp_completed.propensityScores_IPW$w_att),],
    weights = ~ w_att
    )
  
## Construct a table
tabWeighted.w_att <- 
  svyCreateTableOne(
    vars = vars,
    strata = "C3_2w_71", 
    data = res.svydesign.w_att, 
    test = FALSE
    )
## Show table with SMD

sink(
  "results/TableOne.w_att_weighted.txt"
  )
print(tabWeighted.w_att, smd = TRUE)
sink()


## Construct a data frame containing variable name and SMD from all methods
dataPlot <- data.frame(
  variable  = rownames(ExtractSmd(tabUnmatched)),
  rawdata = as.numeric(ExtractSmd(tabUnmatched)),
  weighted_data = as.numeric(ExtractSmd(tabWeighted.w_att))
  )

## Create long-format data for ggplot2
dataPlotMelt <-
  melt(
    data          = dataPlot,
    id.vars       = c("variable"),
    variable.name = "Method",
    value.name    = "SMD"
    )

## Order variable names by magnitude of SMD
varNames <- unique(
  as.character(dataPlot$variable)[
    order(dataPlot$rawdata)
    ]
  )


## Order factor levels in the same order
dataPlotMelt$variable <- 
  factor(
    dataPlotMelt$variable,
    levels = varNames
    )

## Plot using ggplot2

quartz(
  family = "Arial",type = "pdf",
  file =   "results/smd.pdf")
ggplot(
  data = dataPlotMelt[dataPlotMelt$variable %in% c("infection14", "IVIg", "PCT_pre"),],
  mapping = aes(x = variable, y = SMD, group = Method, color = Method)
  ) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
  coord_flip() +
  theme_bw() + theme(legend.key = element_blank())
ggplot(
  data = dataPlotMelt,
  mapping = aes(x = variable, y = SMD, group = Method, color = Method)
  ) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
  coord_flip() +
  theme_bw() + theme(legend.key = element_blank())
dev.off()

# Estimation of the ATT using IPW -------------
# 
# treatment_effect(
#   LDLT_Comp_completed.propensityScores_IPW$C3_2w_71,
#   LDLT_Comp_completed.propensityScores_IPW$mortality90,
#   LDLT_Comp_completed.propensityScores_IPW$w_atc
#   )


## Weighted model (use svyglm function) and without-weight model (use glm function).

#' From help(svyglm):
#' Note
#' svyglm always returns 'model-robust' standard errors;
#'  the Horvitz-Thompson-type standard errors used 
#'  everywhere in the survey package are a generalisation
#'   of the model-robust 'sandwich' estimators. 
#'  In particular, a quasi-Poisson svyglm will return
#'   correct standard errors for relative risk regression models.

res.glmWeighted <-
  svyglm(
    formula = (mortality90 == 1) ~ C3_2w_71,
    family  = quasibinomial(link = "logit"),                  
    design    = res.svydesign.w_att
    )

res.glmOrd <-
  glm(
    formula = (mortality90 == 1) ~ C3_2w_71,
    family  = binomial(link = "logit"),                  
    #design    = res.svydesign.w_att
    data=LDLT_Comp_completed.propensityScores_IPW
  )

res.glmOrd_weighted <-
  glm(
    formula = (mortality90 == 1) ~ C3_2w_71,
    family  = quasibinomial(link = "logit"),                  
    #design    = res.svydesign.w_att
    data=LDLT_Comp_completed.propensityScores_IPW,
    weights = w_att
  )


summary.res.glmWeighted <-
  summary(res.glmWeighted)

print(exp(summary.res.glmWeighted$coefficients[,"Estimate"]))
print(exp(confint(res.glmWeighted)))

#```