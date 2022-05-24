#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Loading Packages
library(shiny)
library(tidyverse)
library(ggsci)
library(shinysky)
library(shinycssloaders)
library(Hmisc)
library(extrafont)
library(dplyr)
library(feather)
loadfonts(device = "pdf")
#Loading and manipulating data
dat <- read_feather(path = "Rlog_for_RNA-Seq-App.feather")
# Getting a vector of gene names to populate autocomplete
colnames(dat)[1] <- "gene"
genes = unique(data.frame(value = dat$gene, label = dat$gene))

# Defining function to calculate SEM
sem <- function(x) sd(x)/sqrt(length(x))
get_geneinfo <- function(x) filter(dat, gene %in% x) %>% top_n(n = 1, wt = gene)

# Define UI for application that plots
ui <- fluidPage(
  
  # Application title
  titlePanel("Gene Expression in the Chicken Neural Crest Over Time"),
  
  # Sidebar menu
  sidebarLayout(
    sidebarPanel(
      # Using selectize input to get list with autocomplete for all valid options (genes). Choices
      # has to be NULL here because the list is huge and would crash the html page if we tried to load it
      selectizeInput(inputId = "g",label="Gene", choices = NULL, multiple = TRUE),
      checkboxInput(inputId = "c", label = "Display Whole Embryo Counts", value = FALSE),
      checkboxInput(inputId = "f", label = "Facet by Gene", value = FALSE),
      checkboxInput(inputId = "ci", label = "Loess Smoothed 80% Confidence Interval", value = FALSE),
      checkboxInput(inputId = "p", label = "Display data points", value = FALSE),
      checkboxInput(inputId = "eb", label = "Display error bars", value = FALSE),
      checkboxInput(inputId = "log", label = "Display Log2FC Enrichment in NC", value = TRUE),
      downloadLink("save", "Download Plot")
    ),
    
    # Plot (see below)
    mainPanel(
      plotOutput("expressionPlot") %>% withSpinner()
    )
  )
)

# Define server logic for plotting
server <- function(input, output, session) {
  options(shiny.usecairo=T)
  updateSelectizeInput(session, inputId = 'g', choices = genes, server = TRUE, options = list(create=FALSE))
  
  output$expressionPlot <- renderPlot(height = 350, {
    
    # Plot the expression plot
    g <- ggplot(
      if (!input$c) {
        filter(dat, condition != 'WE' & gene %in% input$g)
      } else {
        filter(dat, gene %in% input$g)
      },
      aes(x = time,
          y = Rlog,
          color = gene,
          fill = gene
      )
    )
    g <- g + 
      
      {if (length(input$g) != 0) stat_summary(fun = mean, geom = 'line', size=1.5)} +
      
      
      {if (input$log & length(input$g) == 1)
        annotate("text", label = paste0("Log2FC ", round(get_geneinfo(input$g[[1]])[1,2],2),
                                        " Padj ", round(get_geneinfo(input$g[[1]])[1,3],4)),
                 x = 11, y = as.numeric(max(get_geneinfo(input$g[[1]])[,7])) + 2) 
        }+
      
      {if (input$c & input$f & length(input$g) != 0) facet_grid(condition ~ gene, scales = 'fixed')}+
      
      {if (input$c & !input$f & length(input$g) != 0) facet_grid(condition ~ ., scales = 'fixed')}+
      
      {if (!input$c & input$f & length(input$g) != 0) facet_grid(. ~ gene, scales = 'fixed')} +
      
      
      
      {if (input$p & length(input$g) != 0) geom_point(size = 2)}+
      
      
      {if (input$ci & length(input$g) != 0) suppressMessages(stat_smooth(aes(color = NA),
                                                                         method = "loess", formula = "y ~ x",
                                                                         geom = "ribbon", level = 0.8, alpha = 0.25))}+
      
      {if (input$eb & length(input$g) != 0) stat_summary(fun.data = mean_se, geom = "errorbar")} +
      
      
      scale_fill_npg()+
      
      scale_color_npg()+
      
      scale_x_continuous(
        limits = c(6,16),
        breaks = c(6,8,10,12,14,16),
        labels = c('HH6', 'HH8', 'HH10', 'HH12', 'HH14', 'HH16')
      )+
      
      xlab(NULL) +
      
      ylab('Rlog_CPM') +
      
      theme(
        text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 14, face = 'bold'),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, size = 2, color = 'black'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 1, color = 'black'),
        axis.text = element_text(size = 12, color = 'black'),
        panel.grid.major.y = element_line(colour = "gray")
        
      )
    ggsave("plot.pdf")
    g
  })
  output$save <- downloadHandler(
    filename = function(){"plot.pdf"},
    content = function(file){
      file.copy("plot.pdf",file, overwrite = T)
    }
  )

}
# Run the application 
shinyApp(ui = ui, server = server)