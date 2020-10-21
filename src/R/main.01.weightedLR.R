#' Weighted logistec regression.
#' 2020/10/15

Bibtex <- TRUE

dir.sub <- "./src/R/sub"
fn.sub <- list.files(dir.sub)

for(fn.sub.i in fn.sub) source(sprintf("%s/%s",dir.sub,fn.sub.i))

##```
  
load(
  file = sprintf("%s/%s", dir.RData, fn.imported_data)
  )

df.working_data <-
  df.imported_data.completed

df.working_data[,
  var.exposure
  ] <-
  df.working_data[,
    var.exposure.conti
    ] <
  var.cutoff[
    var.cutoff$col_name == var.exposure.conti,
    "cutoff"
    ]

df.working_data[, var.outcome] <-
  df.working_data[,
    var.cutoff[
      var.cutoff$col_name == var.outcome.1,
      "col_name"
      ]
    ] <=
  var.cutoff[
    var.cutoff$col_name == var.outcome.1,"cutoff"
    ]


# Propensity score model ----------

fml.ps_model <-
  sprintf(
    'as.factor(%s) ~ %s',
    var.exposure,
    paste(
      var.Psmodel,
      collapse = "+"
      )
    )


# Finishing coreating analysis data ---------------------------------------

df.working_data_completed <-
  df.working_data[
    !is.na(df.working_data[,var.exposure.conti]),
    ]

## Construct a table
tabUnmatched <-
  CreateTableOne(
    vars = var.smd, 
    strata = var.exposure,
    data = df.working_data_completed, 
    test = FALSE
    )
## Show table with SMD
sink(
  "output/TableOne.txt"
  )
print(tabUnmatched, smd = TRUE)
sink()


# Propensity score model ---------------
propensityScoreModel <-
  glm(
    fml.ps_model,
    family  = binomial(link = "logit"),
    data    = df.working_data_completed,
    na.action = na.exclude
    )

df.working_data_completed.propensityScores <-
  df.working_data_completed %>%
  data.frame()

df.working_data_completed.propensityScores$propensity_score <- 
  predict(
    propensityScoreModel,
    type= "response",
    # type= "prob",
    na.action = na.exclude()
  ) %>%
  unlist()

df.working_data_completed.propensityScores_IPW <-
  IPW_weights(
    treatment = 
      as.numeric(
        df.working_data_completed.propensityScores[,var.exposure]
        ),
    propensity_score = 
      df.working_data_completed.propensityScores$propensity_score,
    dat = df.working_data_completed.propensityScores
    ) %>%
  dplyr::filter(
    !is.na(propensity_score)
    )

res.roc.propensity_score <- roc(
  response = 
    as.factor(
      df.working_data_completed.propensityScores[,var.outcome]
      ),
  predictor = 
    df.working_data_completed.propensityScores$propensity_score
  )

# Plot distribution of propensity score and weighted counts. ----------------
#'
#'

df.working_data_completed.propensityScores_IPW$y_for_dots <-
  (
    as.numeric(
      df.working_data_completed.propensityScores_IPW[,var.exposure]
      )
    )

df.working_data_completed.propensityScores_IPW$var.exposure <-
  df.working_data_completed.propensityScores_IPW[,var.exposure]
df.working_data_completed.propensityScores_IPW$factor.var.exposure <-
  factor(
    df.working_data_completed.propensityScores_IPW[,var.exposure],
    levels = c(FALSE,TRUE),
    unlist(
      strsplit(
        gsub("^\\{(.+)\\}$","\\1",var.label[var.label$col_name==var.exposure,"var.label"]),
        split = "\\}\\{"
        )
      )
    )
  

ggdata.propensityScores <-
  ggplot(
    data =
      df.working_data_completed.propensityScores_IPW,
    aes(
      x = propensity_score
      )
    )

ggdata.propensityScores.weighted_count <-
  ggplot(
    data =
      df.working_data_completed.propensityScores_IPW %>%
      pivot_longer(
        cols =
          c(starts_with("w_at")),
        values_to = "weight",
        names_to = "target_pop"
      ) %>%
      dplyr::filter(
        target_pop == "w_att"
        ),
    aes(
      x = propensity_score,
      weight = weight
    )
  )

pdf(
  file = 
    "output/IPWcount.pdf",
  width = 21
  )
plot(
  ggdata.propensityScores + 
    geom_density(
      aes(
        fill = factor.var.exposure
        ),
      bw="SJ",
#      binwidth = FD,
      alpha=0.5,
      position="identity"
    ) +
  geom_point(
    aes(
      y=y_for_dots,
      x=propensity_score,
      color=factor.var.exposure
      ),
    size=3
  ) +
  theme_bw()
)

ggdata.propensityScores + 
  stat_ecdf(
    aes(
      color = factor.var.exposure
      ),
    alpha=0.5,
    position="identity"
    ) +
  geom_point(
    aes(
      y=y_for_dots,
      x=propensity_score,
      color=factor.var.exposure
      ),
    size=3
    ) +
  theme_bw()

plot(
  ggdata.propensityScores.weighted_count +
    geom_density(
      aes(
        fill= factor.var.exposure
        ),
      bw="SJ",
      alpha=0.5
      ) +
    facet_grid(~target_pop) +
    theme_bw()
  )
plot(
  res.roc.propensity_score,
  print.thres=TRUE 
  )
legend(
  x = 0.6, y=0.5,
  cex = 0.7, 
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
    data =
      df.working_data_completed.propensityScores_IPW[
        !is.na(df.working_data_completed.propensityScores_IPW$propensity_score),
        ],
    weights = ~ w_att
    )
  
## Construct a table
tabWeighted.w_att <- 
  svyCreateTableOne(
    vars = var.smd,
    strata = var.exposure, 
    data = res.svydesign.w_att, 
    test = FALSE
    )
## Show table with SMD

sink(
  "output/TableOne.w_att_weighted.txt"
  )
print(tabWeighted.w_att, smd = TRUE)
sink()


## Construct a data frame containing variable name and SMD from all methods
dataPlot <- 
  data.frame(
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
    ) %>%
  left_join(
    df.col_info,
    by=c("variable"="col_name")
    )

## Order variable names by magnitude of SMD
varNames <- unique(
  as.character(
    dataPlotMelt[
      dataPlotMelt$Method=="rawdata",
      "col_label"]
    )[
    order(
      dataPlotMelt[
        dataPlotMelt$Method=="rawdata",
        "SMD"]
      )
    ]
  )


## Order factor levels in the same order
dataPlotMelt$col_label <- 
  factor(
    dataPlotMelt$col_label,
    levels = varNames
    )

## Plot using ggplot2

quartz(
  family = "Arial",type = "pdf",
  file =   "output/smd.pdf")
ggplot(
  data = dataPlotMelt,
  mapping = aes(x = col_label, y = SMD, group = Method, color = Method)
  ) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
  coord_flip() +
  theme_bw() + theme(legend.key = element_blank())
ggplot(
  data = dataPlotMelt[dataPlotMelt$variable %in% var.Psmodel,],
  mapping = aes(x = col_label, y = SMD, group = Method, color = Method)
  ) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
  coord_flip() +
  theme_bw() + theme(legend.key = element_blank())
dev.off()


## Weighted model (use svyglm function) and without-weight model (use glm function).

#' From help(svyglm):
#' Note
#' svyglm always returns 'model-robust' standard errors;
#'  the Horvitz-Thompson-type standard errors used 
#'  everywhere in the survey package are a generalisation
#'   of the model-robust 'sandwich' estimators. 
#'  In particular, a quasi-Poisson svyglm will return
#'   correct standard errors for relative risk regression models.

fml.gml <- sprintf(
  "as.factor(%s) ~ %s",
  var.outcome,
  var.exposure
  )

res.glmWeighted <-
  svyglm(
    formula = fml.gml,
    family  = quasibinomial(link = "logit"),                  
    design  = res.svydesign.w_att
    )

summary.res.glmWeighted <-
  summary(res.glmWeighted)

sink("output/summary.res.glmWeighted.txt")

print(summary.res.glmWeighted)
print(exp(summary.res.glmWeighted$coefficients[,"Estimate"]))
print(exp(confint(res.glmWeighted)))
sink()
#```