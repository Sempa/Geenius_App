library(tidyverse)
library(plyr)
library(ggplot2)
library(data.table)
library(pracma)
library(cubature)
library(Rmisc)
library(gridExtra)
library(RColorBrewer)
library(splines)
library(minpack.lm)
library(shinydashboard)
library(fontawesome)
library(Cairo)
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinyLP)
library(stringi)
library(coda)
library(jsonlite)

source("pct_curves_uniroot-generic.R")
source("likelihood_function - generic.R")
pct_vs_t_logit_cubic_Geenius <- readRDS("data/pct_data_logitcubic_Geenius.rds")
# likelihood data
pr_t_logit_cubic_Geenius <- readRDS("data/pr_t_logit_evaluations_Geenius.rds")
pr_t_loglog_cubic_Geenius <- readRDS("data/pr_t_loglog_evaluations_Geenius.rds")

# pct_plot <- ggplot(
#   data = pct_vs_t_logit_cubic_Geenius, # subset(data_merge, model_cat == model_type[i]),
#   aes(x = pct_vs_t, y = threshold, group = percentile, color = percentile)
# ) +
#   geom_smooth(se = F, span = .95, size = 2.5) +
#   xlab("Time since infection (days)") +
#   ylab("ODn") +
#   theme_bw() +
#   scale_x_continuous(limits = c(0, 610), breaks = c(seq(0, 610, 200)), expand = c(0, 0)) +
#   scale_colour_brewer(palette = "Paired") +
#   coord_cartesian(ylim = c(0, 3 + .02)) +
#   theme(
#     text = element_text(size = 20),
#     plot.title = element_text(hjust = 0.5),
#     axis.line = element_line(colour = "black"), # panel.grid.major = element_blank(),
#     # panel.grid.minor = element_blank(),
#     axis.text = element_text(size = 18),
#     axis.title = element_text(size = 18),
#     panel.background = element_blank(),
#     panel.border = element_blank(),
#     aspect.ratio = 1,
#     plot.margin = unit(c(0, 0, 0, 0), "null")
#   )
# pct_plot

percentile_func <- function(vec_x, pct) {
  x <- sort(!is.na(vec_x), decreasing = FALSE)
  i_th_observation <- round(pct * (length(vec_x) + 1), 0)
  return(x[i_th_observation])
}

# percentile_func(vec_x = c(1, 3, 5, 6, 9, 11, 12, 13, 19, 21, 22, 32, 35, 36, 45, 44, 
#                      55, 68, 79, 80, 81, 88, 90, 91, 92, 100, 112, 113, 114, 
#                      120, 121, 132, 145, 146, 149, 150, 155, 180, 189, 190),
#            pct = 0.2) = 20
# percentile_func(vec_x = c(1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5, 6, 6), pct = .25) = 2

dat <- data.frame(s_id = 1:6, assay_value_th = c(0.1, 0.2, .4, .6, .8, 1), interval_length = c(100, 120, 180, 130, 200, 205))# assay_value_th <- c(.2, .3 ,.5, 1.5, 2, 2.5, 4)
# source("likelihood_function - geenius.R")
complete_dataset <- data.frame(s_id = NA, l = NA, time_t = NA, bigL = NA, assay_value  = NA, int_length = NA)
# assay_value_th_genious <- assay_value_th_genious
GV_interval_step_genious <- 1
for (j in 1:length(dat$s_id)) {
  assay_value_th_genious <- dat$assay_value_th[j]
  GV_interval_length_genious <- dat$interval_length[j]
dat_combine <- data.frame(s_id = NA, l = NA, time_t = NA, bigL = NA, assay_value  = NA, int_length = NA)
# browser()
  # for (i in 1:length(dat$assay_value_th_genious)) {
    # browser()
    likelihood <- likelihood_fun(
      param_datset = likelihood_param_quad_function(
        dat = as.matrix(pr_t_logit_cubic_Geenius),
        target_assay_value = assay_value_th_genious,
        around_assay_value = seq(0.01, 3, .15),
        t_since_ln = seq(0, GV_interval_length_genious, GV_interval_step_genious) # GV_interval_step_genious
      ),
      assay_value = assay_value_th_genious,
      t_since_ln = seq(0, GV_interval_length_genious, GV_interval_step_genious) # GV_interval_step_genious
    ) %>%
      mutate(assay_value = assay_value_th_genious, s_id = dat$s_id[j], int_length = dat$interval_length[j])
    # browser()
    dat_combine <- rbind(dat_combine, likelihood)
  # }
  
  complete_dataset <- rbind(complete_dataset, dat_combine)
}

complete_dataset <- complete_dataset %>%
  filter(!is.na(time_t)) %>%
  arrange(s_id, time_t)

summary_dataset <- complete_dataset %>%
  mutate(id = paste(s_id, assay_value, sep = '_'))%>%
  group_by(id) %>%
  dplyr::summarise(`mode value` = approx(x = bigL, y = time_t, 
                                  xout = max(bigL, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
            `5th percentile` = approx(x = bigL, y = time_t,
                                      xout = percentile_func(bigL, pct = .05), method = "constant", ties = mean, rule = 2)$y,
            `25th percentile` = approx(x = bigL, y = time_t,
                                       xout = percentile_func(bigL, pct = .25), method = "constant", ties = mean, rule = 2)$y,
            `50th percentile` = approx(x = bigL, y = time_t,
                                       xout = percentile_func(bigL, pct = .5), method = "constant", ties = mean, rule = 2)$y,
            `75th percentile` = approx(x = bigL, y = time_t,
                                       xout = percentile_func(bigL, pct = .75), method = "constant", ties = mean, rule = 2)$y,
            `95th percentile` = approx(x = bigL, y = time_t,
                                       xout = percentile_func(bigL, pct = .95), method = "constant", ties = mean, rule = 2)$y,
            `Posterior percentile` = approx(x = bigL, y = time_t,
                                            xout = percentile_func(bigL, pct = .682), method = "constant", ties = mean, rule = 2)$y
            
  )
f_T <- 0.5
cumulative_posterior <- complete_dataset %>%
  dplyr::group_by(s_id) %>%
  mutate(cumsum_posterior = cumsum(bigL)) %>%
  mutate(f_T = f_T,
    t1 = time_t,
    window_size= round(int_length * f_T, 0),
    t1_f_T = ifelse(t1 + window_size <= int_length, t1 + window_size, int_length)
         ) %>%
  mutate(id  = paste(s_id, t1, sep = '_'))

t1FT_records <- cumulative_posterior %>%
  ungroup() %>%
  group_by(id) %>%
  filter(window_size <=time_t) %>%
  ungroup() %>%
  dplyr::select(id, bigL, cumsum_posterior)

merged_dataset <- cumulative_posterior %>%
  mutate(id  = paste(s_id, t1_f_T, sep = '_')) %>%
  right_join(t1FT_records, by = 'id') %>%
  distinct(id, .keep_all = T) %>%
  mutate(window_probs_t1 = cumsum_posterior.y - cumsum_posterior.x) 
summay_dataset_2 <- merged_dataset %>%
  group_by(s_id) %>%
  dplyr::summarise(max_window_prob = approx(x = window_probs_t1, y = time_t, 
                                    xout = max(window_probs_t1, na.rm = T), method = "constant", ties = mean, rule = 2)$y
  )
  
f_p <- .5

plot2 <- ggplot(
  data = merged_dataset,
  aes(x = t1_f_T, y = c_bigL, colour = s_id)
) +
  geom_line(size = 2.0) +
  xlab("Date of infection (days)") +
  ylab("cumulative posterior") +
  theme_bw() +
  scale_x_reverse(expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0)) +
  # scale_discrete_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#FB9A99", "#6A3D9A", "#B2DF8A")) + # , "#6A3D9A"
  coord_fixed(ratio = 1) +
  theme(
    text = element_text(size = 18),
    axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank(),
    panel.border = element_blank(),
    aspect.ratio = 1,
    plot.margin = unit(c(0, 0, 0, 0), "null")
  )

plot2
