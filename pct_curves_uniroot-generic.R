libraries <- c(
  "plyr", "ggplot2", "tidyverse", "data.table", "pracma", "minpack.lm", #"inctools",
  "stringi", "coda", "cubature", "Rmisc", "gridExtra", "RColorBrewer", "splines"#"truncnorm", 
)
for (x in libraries) {
  library(x, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
}

func_logit_squared <- function(t, parameters) {
  1 / (1 + exp(-(parameters[1] + parameters[2] * t + parameters[3] * t^2 # +
  )))
}


func_logit_cubic <- function(t, parameters) {
  parameters[1] - (1 / (1 + exp(-(parameters[2] + parameters[3] * t + parameters[4] * t^2 + parameters[5] * t^3
  ))))
}

func_cloglog_cubic <- function(t, parameters) {
  parameters[1] - (1 - exp(-exp(parameters[2] + parameters[3] * t + parameters[4] * t^2 + parameters[5] * t^3)))
}

pct_vs_t_fun_logit_quad <- function(data_set, assay_valueTh, time_t, percentile) {
  pct_vs_t <- c(0)
  dat <- data_set %>%
    mutate(eddi_1 = eddi, eddi_2 = eddi^2, eddi_3 = eddi^3)
  counter <- 0
  for (i in 1:length(assay_valueTh)) {
    counter <- counter + 1
    dat <- dat %>%
      mutate(recency = ifelse(assay_value <= assay_valueTh[i] & viral_load > 1000 & eddi <= 1000, 1, 0))
    
    model <- suppressWarnings(glm2::glm2(recency ~ 1 + eddi_1 + eddi_2, # + eddi_3
                                         family = stats::binomial(link = "cloglog"),
                                         data = dat,
                                         control = stats::glm.control(
                                           epsilon = 1e-08,
                                           maxit = 5e4,
                                           trace = FALSE
                                         )
    ))
    for (i in 1:length(percentile)) {
      solve_4_t <- (uniroot(
        f = func_logit_squared, parameters = c(
          percentile[i],
          model$coefficients[[1]], model$coefficients[[2]],
          model$coefficients[[3]]#, model$coefficients[[4]]
        ),
        interval = c(0, 1000), extendInt = "yes" # lower = 0, upper = 1000
      ))$root
      pct_vs_t <- c(pct_vs_t, solve_4_t)
    }
  }
  return(pct_vs_t[-1])
}


pct_vs_t_fun_cloglog_cubic <- function(data_set, assay_valueTh, time_t, percentile) {
  pct_vs_t <- c(0)
  dat <- data_set %>%
    mutate(eddi_1 = eddi, eddi_2 = eddi^2, eddi_3 = eddi^3)
  counter <- 0
  for (i in 1:length(assay_valueTh)) {
    counter <- counter + 1
    dat <- dat %>%
      mutate(recency = ifelse(assay_value <= assay_valueTh[i] & viral_load > 1000 & eddi <= 1000, 1, 0))
    
    model <- suppressWarnings(glm2::glm2(recency ~ 1 + eddi_1 + eddi_2 + eddi_3, #
                                         family = stats::binomial(link = "cloglog"),
                                         data = dat,
                                         control = stats::glm.control(
                                           epsilon = 1e-08,
                                           maxit = 5e4,
                                           trace = FALSE
                                         )
    ))
    for (i in 1:length(percentile)) {
      solve_4_t <- (uniroot(
        f = func_cloglog_cubic, parameters = c(
          percentile[i],
          model$coefficients[[1]], model$coefficients[[2]],
          model$coefficients[[3]], model$coefficients[[4]]
        ),
        interval = c(0, 1000), extendInt = "yes" # lower = 0, upper = 1000
      ))$root
      pct_vs_t <- c(pct_vs_t, solve_4_t)
    }
  }
  return(pct_vs_t[-1])
}


pct_vs_t_fun_logit_cubic <- function(data_set, assay_valueTh, time_t, percentile) {
  pct_vs_t <- c(0)
  dat <- data_set %>%
    mutate(eddi_1 = eddi, eddi_2 = eddi^2, eddi_3 = eddi^3)
  counter <- 0
  for (i in 1:length(assay_valueTh)) {
    counter <- counter + 1
    dat <- dat %>%
      mutate(recency = ifelse(assay_value <= assay_valueTh[i] & viral_load > 1000 & eddi <= 1000, 1, 0))
    
    model <- suppressWarnings(glm2::glm2(recency ~ 1 + eddi_1 + eddi_2 + eddi_3, #
                                         family = stats::binomial(link = "logit"),
                                         data = dat
    ))
    for (i in 1:length(percentile)) {
      solve_4_t <- (uniroot(
        f = func_logit_cubic, parameters = c(
          percentile[i],
          model$coefficients[[1]], model$coefficients[[2]],
          model$coefficients[[3]], model$coefficients[[4]]
        ),
        interval = c(0, 1010), extendInt = "yes"
      ))$root
      pct_vs_t <- c(pct_vs_t, solve_4_t)
    }
  }
  return(pct_vs_t[-1])
}
