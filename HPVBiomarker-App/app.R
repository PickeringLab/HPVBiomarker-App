
library(shiny)
library(readxl)
library(tidyverse)
library(sjmisc)
library(shinythemes)
library(DT)
library(ggplot2)
library(gridExtra)

# Define UI for application that draws a histogram
ui <- fluidPage(
    theme = shinytheme("flatly"),

    navbarPage(
        title = "HPV-Biomarker",
        tabPanel("Import data",
                 sidebarLayout(
                     sidebarPanel(
                         width = 4,
                         #hr(),
                         h3("Import Nanostring data"),
                         br(),
                         fileInput("datatab", "Upload Nanostring output (tab separated txt file) ...", 
                                   accept = ".txt")
                     ),
                     mainPanel(
                         dataTableOutput("datatable")
                     )
                 )),
        tabPanel("Filter data",
                 sidebarLayout(
                     sidebarPanel(
                         width = 4,
                         h3("Filter Nanostring data"),
                         br(),
                         p("The analysis only requires the", strong('Probe Name'), " column and the sample columns."),
                         p("Using the drop-down menu below, please select all the columns that are not required in the analysis, so that they can be filtered out."),
                         p("If by mistake you choose a wrong column, click over it on the selection list and press 'Delete'."),
                         selectInput("filter_columns", "Select the columns to be filtered:", "", selected = NULL, multiple = TRUE),
                         #br(),
                         p(em("Note: Please make sure to filter your data properly before generating the results. The results are wrong until this step is properly carried out.")),
                         p("When you have made your selections properly, click the button below to generate the results."),
                         fluidRow(
                             column(12,
                                    actionButton("gen_results", "Generate results", width = "100%", icon = icon("chart-pie"),
                                                 style = "color: white; background-color: #2c3e50; border-color: #2c3e50"))
                         )
                     ),
                     mainPanel(
                         dataTableOutput("filteredtable")
                     )
                 )),
        tabPanel("Results",
                 navlistPanel(widths = c(3, 9),
                              tabPanel("Table",
                                       br(),
                                       fluidRow(
                                           column(11, offset = 1,
                                                  dataTableOutput("resultstable"))
                                       ),
                                       br(),
                                       em(h4("Note: Please make sure to filter your data before extracting the results. The results are wrong until that step is properly carried out."))
                                       
                              ),
                              tabPanel("Plots",
                                       plotOutput("flowchart"),
                                       fluidRow(
                                           column(6,
                                                  plotOutput("scatter1")),
                                           column(6,
                                                  plotOutput("scatter2"))
                                       )
                              )))
        
                 
        ))





# Define server logic required to draw a histogram
server <- function(input, output) {
    
    datainput <- reactive({
        req(input$datatab)
        
        ext <- tools::file_ext(input$datatab$name)
        switch(ext,
               txt = read_delim(input$datatab$datapath),
               validate("Invalid file; Please upload a .xlsx, .csv or .tsv file")
        )
    })
    
    output$datatable <- renderDataTable(datainput(), options = list(pageLength = 10, scrollX = T))
    
    observeEvent(datainput(), {
        column_choices <- colnames(datainput())
        updateSelectInput(inputId = "filter_columns", choices = column_choices)
    })
    
    filtered <- function(data) {
        data_filtered <- data %>%
            filter(`Probe Name` == "CDA" |
                       `Probe Name` == "C9orf172" |
                       `Probe Name` == "RHOD" |
                       `Probe Name` == "LRMP" |
                       `Probe Name` == "UPP1" |
                       `Probe Name` == "HLF") %>%
            select(-c(input$filter_columns))
    }
    
    results <- function(data) {
        data_filtered <- filtered(data)
        data_transformed <- rotate_df(data_filtered, rn = "Samples", cn = TRUE)
        data_scored <- data_transformed %>%
            mutate(`CDA>AJM1` = case_when(CDA > C9orf172 ~ 1,
                                          CDA <= C9orf172 ~ 0),
                   `RHOD>LRMP` = case_when(RHOD > LRMP ~ 1,
                                           RHOD <= LRMP ~ 0),
                   `UPP1>HLF` = case_when(UPP1 > HLF ~ 1,
                                          UPP1 <= HLF ~ 0))
        data_clustered <- data_scored %>%
            mutate(Classification = case_when((`CDA>AJM1` == 1) & (`RHOD>LRMP` == 1) ~ "C1",
                                              (`CDA>AJM1` == 1) & (`RHOD>LRMP` == 0) ~ "C2",
                                              (`CDA>AJM1` == 0) & (`UPP1>HLF` == 1) ~ "C1",
                                              (`CDA>AJM1` == 0) & (`UPP1>HLF` == 0) ~ "C2"))
        return(data_clustered)
    }
    
    output$filteredtable <- renderDataTable(datatable({
        filtered(datainput())
    },
    extensions = 'Buttons',
    options = list(
        paging = TRUE,
        searching = TRUE,
        fixedColumns = TRUE,
        autoWidth = TRUE,
        ordering = FALSE,
        dom = 'Bliftsp',
        buttons = c('csv', 'excel')
    ),
    class = "display"
    ),
    server = FALSE
    )
    
    output$resultstable <- renderDataTable(datatable({
        p$table
    },
    extensions = 'Buttons',
    options = list(
        paging = TRUE,
        searching = TRUE,
        fixedColumns = TRUE,
        scrollX = TRUE,
        autoWidth = TRUE,
        ordering = FALSE,
        dom = 'Bliftsp',
        buttons = c('csv', 'excel')
    ),
    
    class = "display"
    ),
    server = FALSE)
    
    p <- reactiveValues()
    observeEvent(input$gen_results, {
        
        p$table <- results(datainput())[, -c(2:7)]
        
        p$fl <- tibble(x= 1:100, y= 1:100) %>%
            ggplot(aes(x, y)) +
            scale_x_continuous() +
            scale_y_continuous() +
            theme_test() +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank()) +
            geom_rect(xmin = 30, xmax = 70, ymin = 90, ymax = 100, color = 'steelblue',
                      fill = 'steelblue', size = 0.25, alpha = 0.5) +
            annotate('text', x = 50, y = 95.2, label = paste("CDA > AJM1"), size = 5, colour = "white") +
            geom_rect(xmin = 30, xmax = 45, ymin = 80, ymax = 88, color = 'lightsalmon3',
                      fill = 'lightsalmon3', size = 0.25) +
            annotate('text', x = 38, y = 84.2, label = paste("No =", nrow(subset(results(datainput()), `CDA>AJM1` == 0))), size = 4, color = "#2C3E50") +
            geom_rect(xmin = 55, xmax = 70, ymin = 80, ymax = 88, color = 'seagreen3',
                      fill = 'seagreen3', size = 0.25) +
            annotate('text', x = 63, y = 84.2, label = paste("Yes =", nrow(subset(results(datainput()), `CDA>AJM1` == 1))), size = 4, color = "#2C3E50") +
            geom_segment(x = 37.5, xend = 37.5, y = 79.5, yend = 70, size = 0.5, linejoin = "round", lineend = "round", color = "#2C3E50") +
            geom_segment(x = 62.5, xend = 62.5, y = 79.5, yend = 70, size = 0.5, linejoin = "round", lineend = "round", color = "#2C3E50") +
            geom_segment(x = 37.5, xend = 30, y = 70, yend = 70, size = 0.5, linejoin = "round", lineend = "round", color = "#2C3E50") +
            geom_segment(x = 62.5, xend = 70, y = 70, yend = 70, size = 0.5, linejoin = "round", lineend = "round", color = "#2C3E50") +
            geom_segment(x = 30, xend = 30, y = 70, yend = 50.5, size = 0.5, linejoin = "round", lineend = "round",
                         arrow = arrow(length = unit(2, "mm"), type= "closed"), color = "#2C3E50") +
            geom_segment(x = 70, xend = 70, y = 70, yend = 50.5, size = 0.5, linejoin = "round", lineend = "round",
                         arrow = arrow(length = unit(2, "mm"), type= "closed"), color = "#2C3E50") +
            geom_rect(xmin = 15, xmax = 45, ymin = 40, ymax = 50, color = 'steelblue',
                      fill = 'steelblue', size = 0.25) +
            annotate('text', x = 30, y = 45, label = paste("UPP1 > HLF"), size = 5, color = "white") +
            geom_rect(xmin = 55, xmax = 85, ymin = 40, ymax = 50, color = 'steelblue',
                      fill = 'steelblue', size = 0.25) +
            annotate('text', x = 70, y = 45, label = paste("RHOD > LRMP"), size = 5, color = "white") +
            geom_rect(xmin = 15, xmax = 25, ymin = 30, ymax = 38, color = 'lightsalmon3',
                      fill = 'lightsalmon3', size = 0.25) +
            annotate('text', x = 20, y = 34, 
                     label = paste("No =", nrow(subset(results(datainput()), `CDA>AJM1` == 0 & `UPP1>HLF` == 0))), size = 4, color = "#2C3E50") +
            geom_rect(xmin = 35, xmax = 45, ymin = 30, ymax = 38, color = 'seagreen3',
                      fill = 'seagreen3', size = 0.25) +
            annotate('text', x = 40, y = 34, 
                     label = paste("Yes =", nrow(subset(results(datainput()), `CDA>AJM1` == 0 & `UPP1>HLF` == 1))), size = 4, color = "#2C3E50") +
            geom_rect(xmin = 55, xmax = 65, ymin = 30, ymax = 38, color = 'lightsalmon3',
                      fill = 'lightsalmon3', size = 0.25) +
            annotate('text', x = 60, y = 34, 
                     label = paste("No =", nrow(subset(results(datainput()), `CDA>AJM1` == 1 & `RHOD>LRMP` == 0))), size = 4, color = "#2C3E50") +
            geom_rect(xmin = 75, xmax = 85, ymin = 30, ymax = 38, color = 'seagreen3',
                      fill = 'seagreen3', size = 0.25) +
            annotate('text', x = 80, y = 34, 
                     label = paste("Yes =", nrow(subset(results(datainput()), `CDA>AJM1` == 1 & `RHOD>LRMP` == 1))), size = 4, color = "#2C3E50") +
            geom_segment(x = 20, xend = 20, y = 29.5, yend = 20, size = 0.5, linejoin = "round", lineend = "round",
                         arrow = arrow(length = unit(2, "mm"), type= "closed"), color = "#2C3E50") +
            geom_segment(x = 40, xend = 40, y = 29.5, yend = 20, size = 0.5, linejoin = "round", lineend = "round",
                         arrow = arrow(length = unit(2, "mm"), type= "closed"), color = "#2C3E50") +
            geom_segment(x = 60, xend = 60, y = 29.5, yend = 20, size = 0.5, linejoin = "round", lineend = "round",
                         arrow = arrow(length = unit(2, "mm"), type= "closed"), color = "#2C3E50") +
            geom_segment(x = 80, xend = 80, y = 29.5, yend = 20, size = 0.5, linejoin = "round", lineend = "round",
                         arrow = arrow(length = unit(2, "mm"), type= "closed"), color = "#2C3E50") +
            annotate('text', x = 20, y = 15, label = "C2", size = 5, color = "turquoise3") +
            annotate('text', x = 40, y = 15, label = "C1", size = 5, color = "coral2") +
            annotate('text', x = 60, y = 15, label = "C2", size = 5, color = "turquoise3") +
            annotate('text', x = 80, y = 15, label = "C1", size = 5, color = "coral2") +
            geom_rect(xmin = 30, xmax = 70, ymin = 0, ymax = 10, color = 'steelblue',
                      fill = 'steelblue', size = 0.25) +
            annotate('text', x = 50, y = 5, label = paste("Total samples = ", nrow(results(datainput()))), size = 5, color = "white") +
            labs(title = "DECISION TREE") +
            theme(panel.border = element_blank(),
                  panel.background = element_rect(fill = "white", colour = "white"),
                  plot.title = element_text(face = "bold", hjust = 0.5, color = "#2C3E50"))
        
        p$scatter1 <- ggplot(subset(results(datainput()), `CDA>AJM1` == 0), aes(x = UPP1, y = HLF, color = Classification)) +
            geom_point() +
            geom_abline(color = "#2C3E50") +
            theme_bw() +
            labs(title = "CDA > AJM1") +
            theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#2C3E50"), 
                  axis.title = element_text(face = "bold", color = "#2C3E50"), 
                  axis.text = element_text(color = "#2C3E50"),
                  legend.title = element_blank(),
                  legend.text = element_text(color = "#2C3E50"),
                  panel.background = element_rect(fill = "white", colour = "#2C3E50"),
                  panel.border = element_blank(),
                  plot.background = element_rect(fill = "white"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
        
        p$scatter2 <- ggplot(subset(results(datainput()), `CDA>AJM1` == 1), aes(x = RHOD, y = LRMP, color = Classification)) +
            geom_point() +
            geom_abline(color = "#2C3E50") +
            theme_bw() +
            labs(title = "CDA <= AJM1") +
            theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#2C3E50"), 
                  axis.title = element_text(face = "bold", color = "#2C3E50"), 
                  axis.text = element_text(color = "#2C3E50"),
                  legend.title = element_blank(),
                  legend.text = element_text(color = "#2C3E50"),
                  panel.background = element_rect(fill = "white", colour = "#2C3E50"),
                  panel.border = element_blank(),
                  plot.background = element_rect(fill = "white"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
    })
    
    output$flowchart <- renderPlot({
        p$fl
    }, res = 96)
    
    output$scatter1 <- renderPlot({
        p$scatter1
    }, res = 96)
    
    output$scatter2 <- renderPlot({
        p$scatter2
    }, res = 96)

    
}

# Run the application 
shinyApp(ui = ui, server = server)
