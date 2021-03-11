#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(Seurat)

rdsfiles <- list.files(path = "./data/", pattern = "\\.rds$")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel( title = h1("Feature plot from scRNA seq data"),
                        h3("From pre-processed Seurat object")
               ),

    # Sidebar with a select input 
    sidebarLayout(
        sidebarPanel(
            # Radio button
            radioButtons(inputId = 'choice',
                         label = 'Do you want to precess pre installed data or your own data ?',
                         choices = c("Pre-installed", "Own data"),
                         selected = "Pre-installed"
            ),
            
            conditionalPanel( condition = "input.choice == 'Pre-installed'",
                              # Select pre-installed dataset
                              selectInput(inputId = "piDS",
                                          label = "Dataset",
                                          choices = rdsfiles)
            ),
            
            conditionalPanel( condition = "input.choice == 'Own data'",
                              # Select dataset as .Rds
                              fileInput(inputId = "ownDS",
                                        label = "Dataset",
                                        accept = c('.rds','.Rds') )
                                  
            ),
            
            # Load action button
            actionButton(inputId = "load",
                         label = "Load"),
            
            selectizeInput(inputId = "genes",
                    label = "Gene:",
                    choices = "",
                    multiple = TRUE
            ),
            
            selectInput(inputId = "red",
                        label = "Projections:",
                        choices = ""
            ),
            
            
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            splitLayout(cellWidths = c("50%","50%"), uiOutput('out_dim'), uiOutput('out_ridge'))
            )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    options(shiny.maxRequestSize=100*1024^2)
    
    datasetInput <- reactive({
        if ( input$choice == "Pre-installed") {
            df <- readRDS(paste0("./data/", input$piDS))
            return(df)
        }
        if ( input$choice == "Own data") {
            df <- readRDS(input$ownDS$datapath)
            return(df)
        }
    })
    
    observeEvent(input$load, {
        updateSelectInput(session = session,
                          inputId = "genes",
                          choices = rownames(datasetInput()),
                          selected = rownames(datasetInput())[1])
        updateSelectInput(session = session,
                          inputId = "red",
                          choices = Reductions(datasetInput()),
                          selected = Reductions(datasetInput())[1])
    })
    
    output$out_dim <- renderUI({
        out = list()
        
        if (length(input$genes) == 0){return(NULL)}
        for (i in 1:length(input$genes)){
            out[[i]] <- plotOutput(outputId = paste0("plot_dim",i))
        }
        return(out)
    })
    
    observe({
        for (i in 1:length(input$genes)){
            local({
                ii <- i
                output[[paste0('plot_dim',ii)]] <- renderPlot({
                    return(FeaturePlot(datasetInput(), 
                                       features = input$genes[[ii]],
                                       cols = c("lightgrey","red2"), 
                                       reduction = input$red,
                                       label = TRUE,
                                       combine = FALSE))
                })
            })
        }
    })
    
    output$out_ridge <- renderUI({
        out = list()
        
        if (length(input$genes) == 0){return(NULL)}
        for (i in 1:length(input$genes)){
            out[[i]] <- plotOutput(outputId = paste0("plot_ridge",i))
        }
        return(out)
    })
    
    observe({
        for (i in 1:length(input$genes)){
            local({
                ii <- i
                output[[paste0('plot_ridge',ii)]] <- renderPlot({
                    return(RidgePlot(datasetInput(), 
                                     features = input$genes[[ii]],
                                     combine = FALSE))
                })
            })
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
