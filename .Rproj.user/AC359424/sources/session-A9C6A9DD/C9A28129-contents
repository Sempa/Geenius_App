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

dat <- data.frame(s_id = 1:6, assay_value_th = c(0.1, 0.2, .4, .6, .8, 1), interval_length = c(100, 120, 180, 130, 200, 205))# assay_value_th <- c(.2, .3 ,.5, 1.5, 2, 2.5, 4)
# source("likelihood_function - geenius.R")
complete_dataset <- data.frame(s_id = NA, l = NA, time_t = NA, bigL = NA, assay_value  = NA, int_length = NA)
# assay_value_th_genious <- assay_value_th_genious
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
        t_since_ln = seq(0, GV_interval_length_genious, 2) # GV_interval_step_genious
      ),
      assay_value = assay_value_th_genious,
      t_since_ln = seq(0, GV_interval_length_genious, 2) # GV_interval_step_genious
    ) %>%
      mutate(assay_value = assay_value_th_genious, s_id = dat$s_id[j], int_length = dat$interval_length[j])
    # browser()
    dat_combine <- rbind(dat_combine, likelihood)
  # }
  
  complete_dataset <- rbind(complete_dataset, dat_combine)
}

complete_dataset <- complete_dataset %>%
  filter(!is.na(time_t))

summary_dataset <- complete_dataset %>%
  mutate(id = paste(s_id, assay_value, sep = '_'))%>%
  group_by(id) %>%
  dplyr::summarise(`mode value` = approx(x = bigL, y = time_t, 
                                  xout = max(bigL, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
            `5th percentile` = approx(x = bigL, y = time_t,
                                      xout = quantile(bigL, probs = .05, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
            `25th percentile` = approx(x = bigL, y = time_t,
                                       xout = quantile(bigL, probs = .25, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
            `50th percentile` = approx(x = bigL, y = time_t,
                                       xout = quantile(bigL, probs = .5, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
            `75th percentile` = approx(x = bigL, y = time_t,
                                       xout = quantile(bigL, probs = .75, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
            `95th percentile` = approx(x = bigL, y = time_t,
                                       xout = quantile(bigL, probs = .95, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
            `Posterior percentile` = approx(x = bigL, y = time_t,
                                            xout = quantile(bigL, probs = .682, na.rm = T), method = "constant", ties = mean, rule = 2)$y
            
  )
cumulative_auc <- complete_dataset %>%
  dplyr::group_by(s_id) %>%
  mutate(new_param_F = 0.5,
    t1 = time_t,
    x= int_length * new_param_F,
    t1_new_param_FT = ifelse(t1 + x <= int_length, t1 + x, int_length)
         ) %>%
  mutate(id  = paste(s_id, t1, sep = '_'))

t1FT_records <- cumulative_auc %>%
  ungroup() %>%
  group_by(id) %>%
  filter(x <=time_t) %>%
  ungroup() %>%
  dplyr::select(id, bigL)

merged_dataset <- cumulative_auc %>%
  mutate(id  = paste(s_id, t1_new_param_FT, sep = '_')) %>%
  left_join(t1FT_records, by = 'id') %>%
  distinct(id, .keep_all = T) %>%
  mutate(c_bigL = ((t1_new_param_FT- time_t)/2) * (bigL.x + bigL.y),
         s_id = as.character(s_id))
  

plot2 <- ggplot(
  data = merged_dataset,
  aes(x = t1_new_param_FT, y = c_bigL, colour = s_id)
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
