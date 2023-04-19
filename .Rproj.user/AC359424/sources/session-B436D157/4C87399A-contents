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

# setwd("~/GitHub/biomarker-growth-curves/BM_posteriors_pctcurves")
# source("functions_file.r")
source("pct_curves_uniroot-generic.R")
source("likelihood_function - generic.R")
# shinylogs::track_usage(
#   storage_mode = store_custom(path = "\dropbox\logs")
# )

# Percentile curves
pct_vs_t_logit_cubic_Geenius <- readRDS("data/pct_data_logitcubic_Geenius.rds")
# likelihood data
pr_t_logit_cubic_Geenius <- readRDS("data/pr_t_logit_evaluations_Geenius.rds")
pr_t_loglog_cubic_Geenius <- readRDS("data/pr_t_loglog_evaluations_Geenius.rds")

ui <- navbarPage("Biomarker posteriors and percentile curves",
  theme = shinytheme("flatly"),

  tabPanel(
    title = "Geenius avidity assay",
    sidebarLayout(
      sidebarPanel(
        width = 6,
        box(
          title = "", # Geenius avidity assay
          solidHeader = TRUE,
          numericInput(
            inputId = "assay_value_th_genious",
            label = "value of the recency marker",
            min = 0.001,
            max = 2.5,
            value = 0.5,
            step = 0.001
          ),
          numericInput(
            inputId = "GV_interval_length_genious",
            label = "What is the time interval between last negative and first positive",
            min = 50,
            max = 600,
            value = 400,
            step = 1
          ),
          numericInput(
            inputId =  "GV_percentile",#"GV_interval_step_genious",
            label = "Posterior percentile",
            min = 0.1,
            max = 0.999,
            value = 0.682,
            step = 0.001
          ),
        )
      ),
      fluidRow(
        box("Posterior curve", plotOutput("distPlot_genious")),
        box("Percentile curve", plotOutput("pctPlot_genious")),
        box("Summary table", tableOutput("summary_table_genious"))
      )

    )
  )
)


server <- shinyServer(function(input, output) {

  output$pctPlot_genious <- renderPlot({
    pct_plot <- ggplot(
      data = pct_vs_t_logit_cubic_Geenius, # subset(data_merge, model_cat == model_type[i]),
      aes(x = pct_vs_t, y = threshold, group = percentile, color = percentile)
    ) +
      geom_smooth(se = F, span = .95, size = 2.5) +
      xlab("Time since infection (days)") +
      ylab("ODn") +
      theme_bw() +
      scale_x_continuous(limits = c(0, 610), breaks = c(seq(0, 610, 200)), expand = c(0, 0)) +
      scale_colour_brewer(palette = "Paired") +
      coord_cartesian(ylim = c(0, 3 + .02)) +
      theme(
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black"), # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.background = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1,
        plot.margin = unit(c(0, 0, 0, 0), "null")
      )
    pct_plot
  })

  output$distPlot_genious <- renderPlot({
    assay_value_th_genious <- input$assay_value_th_genious
    dat_combine <- data.frame(l = NA, time_t = NA, bigL = NA, assay_value = NA)
    for (i in 1:length(input$assay_value_th_genious)) {
      # browser()
      likelihood <- likelihood_fun(
        param_datset = likelihood_param_quad_function(
          dat = as.matrix(pr_t_logit_cubic_Geenius),
          target_assay_value = input$assay_value_th_genious[i],
          around_assay_value = seq(0.01, 3, .15),
          t_since_ln = seq(0, input$GV_interval_length_genious, 2)#input$GV_interval_step_genious
        ),
        assay_value = input$assay_value_th_genious[i],
        t_since_ln = seq(0, input$GV_interval_length_genious, 2)#input$GV_interval_step_genious
      ) %>%
        mutate(assay_value = input$assay_value_th_genious[i])
      dat_combine <- rbind(dat_combine, likelihood)
    }

    dat_combine <- dat_combine %>%
      filter(!is.na(time_t)) %>%
      mutate(`assay value` = as.factor(assay_value))
    
    plot2 <- ggplot(
      data = dat_combine,
      aes(x = time_t, y = bigL, colour = `assay value`)
    ) +
      geom_line(size = 2.0) +
      xlab("Date of infection (days)") +
      ylab("Posterior Density") +
      theme_bw() +
      scale_x_reverse(expand = c(0, 0)) +
      scale_y_continuous(breaks = c(0)) +
      scale_colour_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#FB9A99")) + # , "#6A3D9A"
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
  })   
  
  output$summary_table_genious = renderTable({
    assay_value_th_genious <- input$assay_value_th_genious
    dat_combine <- data.frame(l = NA, time_t = NA, bigL = NA, assay_value = NA)
    for (i in 1:length(input$assay_value_th_genious)) {
      # browser()
      likelihood <- likelihood_fun(
        param_datset = likelihood_param_quad_function(
          dat = as.matrix(pr_t_logit_cubic_Geenius),
          target_assay_value = input$assay_value_th_genious[i],
          around_assay_value = seq(0.01, 3, .15),
          t_since_ln = seq(0, input$GV_interval_length_genious, 2)#input$GV_interval_step_genious
        ),
        assay_value = input$assay_value_th_genious[i],
        t_since_ln = seq(0, input$GV_interval_length_genious, 2)#input$GV_interval_step_genious
      ) %>%
        mutate(assay_value = input$assay_value_th_genious[i])
      dat_combine <- rbind(dat_combine, likelihood)
    }
    
    dat_combine <- dat_combine %>%
      filter(!is.na(time_t)) %>%
      mutate(`assay value` = as.factor(assay_value))
    
    summary_table_gi_1 <- dat_combine %>%
      filter(!is.na(assay_value)) %>%
      group_by(assay_value) %>%
      summarise(`mode value` = approx(x = bigL, y = time_t, 
                                    xout = max(bigL, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
                `5th percentile` = approx(x = bigL, y = time_t,
                                          xout = quantile(bigL, probs = .05, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
                `25th percentile` = approx(x = bigL, y = time_t,
                                           xout = quantile(bigL, probs = .25, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
                `50th percentile` = approx(x = bigL, y = time_t,
                                           xout = quantile(bigL, probs = .5, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
                `Posterior percentile` = approx(x = bigL, y = time_t,
                                           xout = quantile(bigL, probs = input$GV_percentile, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
                `75th percentile` = approx(x = bigL, y = time_t,
                                           xout = quantile(bigL, probs = .75, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
                `95th percentile` = approx(x = bigL, y = time_t,
                                           xout = quantile(bigL, probs = .95, na.rm = T), method = "constant", ties = mean, rule = 2)$y
                
      )

  })

})

# shinyApp(ui, server)

# library(rsconnect)
# rsconnect::deployApp('path/to/your/app')
# shinyLP::runExample()
# tmp.enc <- options()$encoding
# options(encoding = "UTF-8")
# rsconnect:: deployApp()
# rsconnect::showLogs()
# library(rsconnect)
# options(rsconnect.http.trace = TRUE, rsconnect.error.trace = TRUE, rsconnect.http.verbose = TRUE)
# rsconnect::setAccountInfo(name='sempa',
#                           token='22479A83F2642CA567F4117568727A0E',
#                           secret='LpxqqXJYT1OD+2I4+74Rt88d0lg/Dq5To11nHrHP')
shinyApp(ui = ui, server = server)
