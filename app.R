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
  # tabPanel("home",
  #   icon = icon("home"),
  #   tags$head(tags$script(HTML('
  #       var fakeClick = function(tabName) {
  #         var dropdownList = document.getElementsByTagName("a");
  #         for (var i = 0; i < dropdownList.length; i++) {
  #           var link = dropdownList[i];
  #           if(link.getAttribute("data-value") == tabName) {
  #             link.click();
  #           };
  #         }
  #       };
  #     '))),
  #   basicPage(
  #     uiOutput(outputId = "text")
  #   ),
  #   jumbotron("Welcome to the continuous diagnostic markers' infection dating tool",
  #     "This App develops biomarker posteriors and percentile curve for the following assays",
  #     button = FALSE
  #     
  #   ),
  #   fluidRow(
  #     column(6, panel_div(
  #       class_type = "primary", panel_title = "Sedia LAg Assay",
  #       box("Sedia LAg Assay", onclick = "fakeClick('Sedia LAg Assay')")
  #     )),
  #     column(6, panel_div(
  #       class_type = "success", panel_title = "Maxim LAg Assay",
  #       box("Maxim LAg Assay", onclick = "fakeClick('Maxim LAg Assay')")
  #     ))
  #   ),
  #   fluidRow( 
  #     column(6, panel_div(
  #       class_type = "info", panel_title = "Architect avidity assay",
  #       box("Architect avidity assay", onclick = "fakeClick('Architect avidity assay')")
  #     )),
  #     column(6, panel_div(
  #       class_type = "danger", panel_title = "Architect Unmodified assay",
  #       box("Architect Unmodified assay", onclick = "fakeClick('Architect Unmodified assay')")
  #     ))
  #   ),
  #   fluidRow( # box("Architect avidity assay", onclick = "fakeClick('Architect avidity assay')")
  #     column(6, panel_div(
  #       class_type = "warning", panel_title = "Bio-rad Geenius avidity assay",
  #       box("Geenius avidity assay", onclick = "fakeClick('Geenius avidity assay')")
  #     ))
  #   ) 
  # ),
  # tabPanel(
  #   "Sedia LAg Assay",
  #   sidebarLayout(
  #     sidebarPanel(
  #       width = 6,
  #       box(
  #         title = "", # Sedia LAg Assay
  #         solidHeader = TRUE,
  #         numericInput(
  #           inputId = "assay_value_th_Sedia",
  #           label = "Value of the recency marker",
  #           min = 0.001,
  #           max = 4,
  #           value = 0.5,
  #           step = 0.001
  #         ),
  #         numericInput(
  #           inputId = "GV_interval_length_Sedia",
  #           label = "What is the time interval between last negative and first positive",
  #           min = 50,
  #           max = 600,
  #           value = 400,
  #           step = 1
  #         ),
  #         numericInput(
  #           inputId = "GV_interval_step_Sedia",
  #           label = "Inter-test interval step size",
  #           min = 1,
  #           max = 7,
  #           value = 2,
  #           step = 1
  #         ),
  #       )
  #     ),
  #     fluidRow(
  #       box("Posterior curve", plotOutput("distPlot_Sedia")),
  #       box("Percentile curve", plotOutput("pctPlot_Sedia"))
  #     )
  #     
  #   )
  # )#,
  # tabPanel(
  #   "Maxim LAg Assay",
  #   sidebarLayout(
  #     sidebarPanel(
  #       width = 6,
  #       box(
  #         title = "",
  #         solidHeader = TRUE,
  #         numericInput(
  #           inputId = "assay_value_th_maxim",
  #           label = "Value of the recency marker",
  #           min = 0.001,
  #           max = 4,
  #           value = 0.5,
  #           step = 0.001
  #         ),
  #         numericInput(
  #           inputId = "GV_interval_length_Maxim",
  #           label = "What is the time interval between last negative and first positive",
  #           min = 50,
  #           max = 600,
  #           value = 400,
  #           step = 1
  #         ),
  #         numericInput(
  #           inputId = "GV_interval_step_Maxim",
  #           label = "Inter-test interval step size",
  #           min = 1,
  #           max = 7,
  #           value = 2,
  #           step = 1
  #         ),
  #       )
  #     ),
  #     fluidRow(
  #       box("Posterior curve", plotOutput("distPlot_Maxim")),
  #       box("Percentile curve", plotOutput("pctPlot_Maxim"))
  #     )
  # 
  #   )
  # ),
  # tabPanel(
  #   "Architect avidity assay",
  #   sidebarLayout(
  #     sidebarPanel(
  #       width = 8,
  #       box(
  #         title = "", # Architect avidity assay
  #         solidHeader = TRUE,
  #         numericInput(
  #           inputId = "assay_value_th_architect",
  #           label = "value of the recency marker",
  #           min = 20,
  #           max = 120,
  #           value = 70,
  #           step = .01
  #         ),
  #         numericInput(
  #           inputId = "GV_interval_length_architect",
  #           label = "What is the time interval between last negative and first positive",
  #           min = 50,
  #           max = 600,
  #           value = 400,
  #           step = 1
  #         ),
  #         numericInput(
  #           inputId = "GV_interval_step_architect",
  #           label = "Inter-test interval step size",
  #           min = 1,
  #           max = 7,
  #           value = 2,
  #           step = 1
  #         ),
  #       )
  #     ),
  #     fluidRow(
  #       box("Posterior curve", plotOutput("distPlot_architect")),
  #       box("Percentile curve", plotOutput("pctPlot_architect"))
  #     )
  #     
  #   )
  # ),
  # tabPanel(
  #   "Architect Unmodified assay",
  #   sidebarLayout(
  #     sidebarPanel(
  #       width = 8,
  #       box(
  #         title = "", # Architect Unmodified assay
  #         solidHeader = TRUE,
  #         numericInput(
  #           inputId = "assay_value_th_ArchitectUnmodified",
  #           label = "value of the recency marker",
  #           min = 0,
  #           max = 500,
  #           value = 70,
  #           step = .01
  #         ),
  #         numericInput(
  #           inputId = "GV_interval_length_ArchitectUnmodified",
  #           label = "What is the time interval between last negative and first positive",
  #           min = 50,
  #           max = 600,
  #           value = 400,
  #           step = 1
  #         ),
  #         numericInput(
  #           inputId = "GV_interval_step_ArchitectUnmodified",
  #           label = "Inter-test interval step size",
  #           min = 1,
  #           max = 7,
  #           value = 2,
  #           step = 1
  #         ),
  #       )
  #     ),
  #     fluidRow(
  #       box("Posterior curve", plotOutput("distPlot_ArchitectUnmodified")),
  #       box("Percentile curve", plotOutput("pctPlot_ArchitectUnmodified"))
  #     )
  #     
  #   )
  # ),
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
            inputId = "GV_interval_step_genious",
            label = "Inter-test interval step size",
            min = 1,
            max = 7,
            value = 2,
            step = 1
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
  # output$text <- renderText({
  #   HTML("<b>","Welcome to the continuous diagnostic markers’ infection dating tool","</b>")
  # })
  # output$pctPlot_Sedia <- renderPlot({
  #   pct_plot <- ggplot(
  #     data = pct_vs_t_logit_cubic_sedia,
  #     aes(x = pct_vs_t, y = threshold, group = percentile, color = percentile)
  #   ) +
  #     geom_smooth(se = F, span = .3, size = 2.5) +
  #     xlab("Time since infection (days)") +
  #     ylab("ODn") +
  #     theme_bw() +
  #     scale_x_continuous(limits = c(0, 610), breaks = c(seq(0, 610, 200)), expand = c(0, 0)) +
  #     scale_colour_brewer(palette = "Paired") +
  #     coord_cartesian(ylim = c(0, 4 + .02)) +
  #     theme(
  #       text = element_text(size = 20),
  #       plot.title = element_text(hjust = 0.5),
  #       axis.line = element_line(colour = "black"),
  #       axis.text = element_text(size = 18),
  #       axis.title = element_text(size = 18),
  #       panel.background = element_blank(),
  #       panel.border = element_blank(),
  #       aspect.ratio = 1,
  #       plot.margin = unit(c(0, 0, 0, 0), "null")
  #     )
  #   pct_plot
  # })
  # 
  # output$distPlot_Sedia <- renderPlot({
  #   assay_value_th_Sedia <- input$assay_value_th_Sedia
  #   dat_combine <- data.frame(l = NA, time_t = NA, bigL = NA, assay_value = NA)
  #   for (i in 1:length(input$assay_value_th_Sedia)) {
  #     likelihood <- likelihood_fun(
  #       param_datset = likelihood_param_quad_function(
  #         dat = as.matrix(pr_t_logit_cubic_sedia),
  #         target_assay_value = input$assay_value_th_Sedia[i],
  #         around_assay_value = seq(0.01, 5, .15),
  #         t_since_ln = seq(0, input$GV_interval_length_Sedia, input$GV_interval_step_Sedia)
  #       ),
  #       assay_value = input$assay_value_th_Sedia[i],
  #       t_since_ln = seq(0, input$GV_interval_length_Sedia, input$GV_interval_step_Sedia)
  #     ) %>%
  #       mutate(assay_value = input$assay_value_th_Sedia[i])
  #     dat_combine <- rbind(dat_combine, likelihood)
  #   }
  # 
  #   dat_combine <- dat_combine %>%
  #     filter(!is.na(time_t)) %>%
  #     mutate(`assay value` = as.factor(assay_value))
  #   plot2 <- ggplot(
  #     data = dat_combine,
  #     aes(x = time_t, y = bigL, colour = `assay value`)
  #   ) +
  #     geom_line(size = 2.0) +
  #     xlab("Date of infection (days)") +
  #     ylab("Posterior Density") +
  #     theme_bw() +
  #     scale_x_reverse(expand = c(0, 0)) +
  #     scale_y_continuous(breaks = c(0)) +
  #     scale_colour_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#FB9A99")) + # , "#6A3D9A"
  #     coord_fixed(ratio = 1) +
  #     theme(
  #       text = element_text(size = 18),
  #       axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(), panel.background = element_blank(),
  #       panel.border = element_blank(),
  #       aspect.ratio = 1,
  #       plot.margin = unit(c(0, 0, 0, 0), "null")
  #     )
  # 
  #   plot2
  # })
  # output$pctPlot_Maxim <- renderPlot({
  #   pct_plot <- ggplot(
  #     data = pct_vs_t_logit_cubic_Maxim,
  #     aes(x = pct_vs_t, y = threshold, group = percentile, color = percentile)
  #   ) +
  #     geom_smooth(se = F, span = .6, size = 2.5) +
  #     xlab("Time since infection (days)") +
  #     ylab("ODn") +
  #     theme_bw() +
  #     scale_x_continuous(limits = c(0, 610), breaks = c(seq(0, 610, 200)), expand = c(0, 0)) +
  #     scale_colour_brewer(palette = "Paired") +
  #     coord_cartesian(ylim = c(0, 3.8 + .02)) +
  #     theme(
  #       text = element_text(size = 20),
  #       plot.title = element_text(hjust = 0.5),
  #       axis.line = element_line(colour = "black"),
  #       # panel.grid.minor = element_blank(),
  #       axis.text = element_text(size = 18),
  #       axis.title = element_text(size = 18),
  #       panel.background = element_blank(),
  #       panel.border = element_blank(),
  #       aspect.ratio = 1,
  #       plot.margin = unit(c(0, 0, 0, 0), "null")
  #     )
  #   pct_plot
  # })
  # 
  # output$distPlot_Maxim <- renderPlot({
  #   assay_value_th_maxim <- input$assay_value_th_maxim
  #   dat_combine <- data.frame(l = NA, time_t = NA, bigL = NA, assay_value = NA)
  #   for (i in 1:length(input$assay_value_th_maxim)) {
  #     likelihood <- likelihood_fun(
  #       param_datset = likelihood_param_quad_function(
  #         dat = as.matrix(pr_t_logit_cubic_Maxim),
  #         target_assay_value = input$assay_value_th_maxim[i],
  #         around_assay_value = seq(0.01, 5, .15),
  #         t_since_ln = seq(0, input$GV_interval_length_Maxim, input$GV_interval_step_Maxim)
  #       ),
  #       assay_value = input$assay_value_th_maxim[i],
  #       t_since_ln = seq(0, input$GV_interval_length_Maxim, input$GV_interval_step_Maxim)
  #     ) %>%
  #       mutate(assay_value = input$assay_value_th_maxim[i])
  #     dat_combine <- rbind(dat_combine, likelihood)
  #   }
  # 
  #   dat_combine <- dat_combine %>%
  #     filter(!is.na(time_t)) %>%
  #     mutate(`assay value` = as.factor(assay_value))
  #   plot2 <- ggplot(
  #     data = dat_combine,
  #     aes(x = time_t, y = bigL, colour = `assay value`)
  #   ) +
  #     geom_line(size = 2.0) +
  #     xlab("Date of infection (days)") +
  #     ylab("Posterior Density") +
  #     theme_bw() +
  #     scale_x_reverse(expand = c(0, 0)) +
  #     scale_y_continuous(breaks = c(0)) +
  #     scale_colour_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#FB9A99")) + # , "#6A3D9A"
  #     coord_fixed(ratio = 1) +
  #     theme(
  #       text = element_text(size = 18),
  #       axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(), panel.background = element_blank(),
  #       panel.border = element_blank(),
  #       aspect.ratio = 1,
  #       plot.margin = unit(c(0, 0, 0, 0), "null")
  #     )
  # 
  #   plot2
  # })
  # output$pctPlot_architect <- renderPlot({
  #   pct_plot <- ggplot(
  #     data = pct_vs_t_logit_cubic_architect,
  #     aes(x = pct_vs_t, y = threshold, group = percentile, color = percentile)
  #   ) +
  #     geom_smooth(se = F, span = .90, size = 2.0) +
  #     xlab("Time since infection (days)") +
  #     ylab("ODn") +
  #     theme_bw() +
  #     scale_x_continuous(limits = c(0, 610), breaks = c(seq(0, 610, 200)), expand = c(0, 0)) +
  #     scale_colour_brewer(palette = "Paired") +
  #     coord_cartesian(ylim = c(20, 140 + .02)) +
  #     theme(
  #       text = element_text(size = 20),
  #       plot.title = element_text(hjust = 0.5),
  #       axis.line = element_line(colour = "black"),
  #       axis.text = element_text(size = 18),
  #       axis.title = element_text(size = 18),
  #       panel.background = element_blank(),
  #       panel.border = element_blank(),
  #       aspect.ratio = 1,
  #       plot.margin = unit(c(0, 0, 0, 0), "null")
  #     )
  #   pct_plot
  # })
  # 
  # output$distPlot_architect <- renderPlot({
  #   assay_value_th_architect <- input$assay_value_th_architect
  #   dat_combine <- data.frame(l = NA, time_t = NA, bigL = NA, assay_value = NA)
  #   for (i in 1:length(input$assay_value_th_architect)) {
  #     likelihood <- likelihood_fun(
  #       param_datset = likelihood_param_quad_function(
  #         dat = as.matrix(pr_t_logit_cubic_architect),
  #         target_assay_value = input$assay_value_th_architect[i],
  #         around_assay_value = seq(20, 245, .15),
  #         t_since_ln = seq(0, input$GV_interval_length_architect, input$GV_interval_step_architect)
  #       ),
  #       assay_value = input$assay_value_th_architect[i],
  #       t_since_ln = seq(0, input$GV_interval_length_architect, input$GV_interval_step_architect)
  #     ) %>%
  #       mutate(assay_value = input$assay_value_th_architect[i])
  #     dat_combine <- rbind(dat_combine, likelihood)
  #   }
  # 
  #   dat_combine <- dat_combine %>%
  #     filter(!is.na(time_t)) %>%
  #     mutate(`assay value` = as.factor(assay_value))
  #   plot2 <- ggplot(
  #     data = dat_combine,
  #     aes(x = time_t, y = bigL, colour = `assay value`)
  #   ) +
  #     geom_line(size = 2.0) +
  #     xlab("Date of infection (days)") +
  #     ylab("Posterior Density") +
  #     theme_bw() +
  #     scale_x_reverse(expand = c(0, 0)) +
  #     scale_y_continuous(breaks = c(0)) +
  #     scale_colour_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#FB9A99")) + # , "#6A3D9A"
  #     coord_fixed(ratio = 1) +
  #     theme(
  #       text = element_text(size = 18),
  #       axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(), panel.background = element_blank(),
  #       panel.border = element_blank(),
  #       aspect.ratio = 1,
  #       plot.margin = unit(c(0, 0, 0, 0), "null")
  #     )
  # 
  #   plot2
  # })
  # output$pctPlot_ArchitectUnmodified <- renderPlot({
  #   pct_plot <- ggplot(
  #     data = pct_vs_t_logit_cubic_ArchitectUnmodified,
  #     aes(x = pct_vs_t, y = threshold, group = percentile, color = percentile)
  #   ) +
  #     geom_smooth(se = F, span = .55, size = 2.0) +
  #     xlab("Time since infection (days)") +
  #     ylab("ODn") +
  #     theme_bw() +
  #     scale_x_continuous(limits = c(0, 610), breaks = c(seq(0, 610, 200)), expand = c(0, 0)) +
  #     scale_colour_brewer(palette = "Paired") +
  #     coord_cartesian(ylim = c(0, 600 + .02)) +
  #     theme(
  #       text = element_text(size = 20),
  #       plot.title = element_text(hjust = 0.5),
  #       axis.line = element_line(colour = "black"),
  #       axis.text = element_text(size = 18),
  #       axis.title = element_text(size = 18),
  #       panel.background = element_blank(),
  #       panel.border = element_blank(),
  #       aspect.ratio = 1,
  #       plot.margin = unit(c(0, 0, 0, 0), "null")
  #     )
  #   pct_plot
  # })
  # 
  # output$distPlot_ArchitectUnmodified <- renderPlot({
  #   assay_value_th_ArchitectUnmodified <- input$assay_value_th_ArchitectUnmodified
  #   dat_combine <- data.frame(l = NA, time_t = NA, bigL = NA, assay_value = NA)
  #   for (i in 1:length(input$assay_value_th_ArchitectUnmodified)) {
  #     likelihood <- likelihood_fun(
  #       param_datset = likelihood_param_quad_function(
  #         dat = as.matrix(pr_t_logit_cubic_ArchitectUnmodified),
  #         target_assay_value = input$assay_value_th_ArchitectUnmodified[i],
  #         around_assay_value = seq(20, 245, 0.15),
  #         t_since_ln = seq(0, input$GV_interval_length_ArchitectUnmodified, input$GV_interval_step_ArchitectUnmodified)
  #       ),
  #       assay_value = input$assay_value_th_ArchitectUnmodified[i],
  #       t_since_ln = seq(0, input$GV_interval_length_ArchitectUnmodified, input$GV_interval_step_ArchitectUnmodified)
  #     ) %>%
  #       mutate(assay_value = input$assay_value_th_ArchitectUnmodified[i])
  #     dat_combine <- rbind(dat_combine, likelihood)
  #   }
  # 
  #   dat_combine <- dat_combine %>%
  #     filter(!is.na(time_t)) %>%
  #     mutate(`assay value` = as.factor(assay_value))
  #   plot2 <- ggplot(
  #     data = dat_combine,
  #     aes(x = time_t, y = bigL, colour = `assay value`)
  #   ) +
  #     geom_line(size = 2.0) +
  #     xlab("Date of infection (days)") +
  #     ylab("Posterior Density") +
  #     theme_bw() +
  #     scale_x_reverse(expand = c(0, 0)) +
  #     scale_y_continuous(breaks = c(0)) +
  #     scale_colour_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#FB9A99")) + # , "#6A3D9A"
  #     coord_fixed(ratio = 1) +
  #     theme(
  #       text = element_text(size = 18),
  #       axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(), panel.background = element_blank(),
  #       panel.border = element_blank(),
  #       aspect.ratio = 1,
  #       plot.margin = unit(c(0, 0, 0, 0), "null")
  #     )
  # 
  #   plot2
  # })
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
          t_since_ln = seq(0, input$GV_interval_length_genious, input$GV_interval_step_genious)
        ),
        assay_value = input$assay_value_th_genious[i],
        t_since_ln = seq(0, input$GV_interval_length_genious, input$GV_interval_step_genious)
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
          t_since_ln = seq(0, input$GV_interval_length_genious, input$GV_interval_step_genious)
        ),
        assay_value = input$assay_value_th_genious[i],
        t_since_ln = seq(0, input$GV_interval_length_genious, input$GV_interval_step_genious)
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
      summarise(`mode val` = approx(x = bigL, y = time_t, 
                                  xout = max(bigL, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
                `25th percentile` = approx(x = bigL, y = time_t,
                                    xout = quantile(bigL, probs = .25, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
                `50th percentile` = approx(x = bigL, y = time_t,
                                    xout = quantile(bigL, probs = .5, na.rm = T), method = "constant", ties = mean, rule = 2)$y,
                `75th percentile` = approx(x = bigL, y = time_t,
                xout = quantile(bigL, probs = .75, na.rm = T), method = "constant", ties = mean, rule = 2)$y
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
