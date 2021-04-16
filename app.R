#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

if (!require(shiny)) { install.packages("shiny") }
if (!require(shinythemes)) { install.packages("shinythemes") }
if (!require(Seurat)) { install.packages("Seurat") }

rdsfiles <- list.files(path = "./data/", pattern = "\\.rds$" , ignore.case = T)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel( title = h1("Feature plot from scRNA seq data"),
                h3("From pre-processed Seurat object")
    ),
    # Theme to use
    themeSelector(),
    
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
                                        accept = c('.rds','.Rds'))
            ),
            
            # Load action button
            actionButton(inputId = "load",
                         label = "Load"),
            
            selectizeInput(inputId = "genes",
                           label = "Features:",
                           choices = "",
                           multiple = TRUE
            ),
            
            selectInput(inputId = "red",
                        label = "Projections:",
                        choices = ""
            ),
            
            selectInput(inputId = "type",
                        label = "Type of features plot:",
                        choices = c("RidgePlot", "VlnPlot", "DotPlot", "HeatMap"),
                        multiple = FALSE
            ),
            
            
            conditionalPanel( condition = "input.type == 'VlnPlot'",
                              # Choose to split or not
                              selectInput(inputId = "split",
                                          label = "Split data by:",
                                          choices = c("Origins", "None"),
                                          multiple = FALSE),
                              conditionalPanel( condition = "input.split != 'None'",
                                                # Choose to split plot or nor
                                                radioButtons(inputId = "split.plot",
                                                             label = "Do you want to split your Violin Plots",
                                                             choices = c("No", "Yes"),
                                                             selected = "No"))
            )
            
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(type = "pills",
                        tabPanel("Unique Feature", splitLayout(cellWidths = c("50%","50%"), uiOutput('out_dim'), uiOutput('out_feat')) ),
                        tabPanel("Multiple Features", uiOutput('out_multi'), uiOutput('out_multi2')),
                        tabPanel("Tables", uiOutput('out_table'), downloadButton('export_table', 'Export csv'))
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    options(shiny.maxRequestSize=1024*1024^2)
    
    datasetInput <- reactive({
        if ( input$choice == "Pre-installed") {
            data.set <- readRDS(paste0("./data/", input$piDS))
            return(data.set)
        }
        if ( input$choice == "Own data") {
            data.set <- readRDS(input$ownDS$datapath)
            return(data.set)
        }
    })
    
    Avg <- reactive({ AverageExpression(datasetInput() , features = input$genes)$RNA})
    
    # AvgMtx <- AverageExpression(datasetInput())
    
    observeEvent(input$load, {
        updateSelectInput(session = session,
                          inputId = "genes",
                          choices = c( rownames(datasetInput()) , colnames(datasetInput()@meta.data) ) ,
                          selected = rownames(datasetInput())[1])
        updateSelectInput(session = session,
                          inputId = "red",
                          choices = Reductions(datasetInput()),
                          selected = Reductions(datasetInput())[1])
        updateSelectInput(session = session,
                          inputId = "split",
                          choices = c("None", colnames(datasetInput()@meta.data)[grep(colnames(datasetInput()@meta.data), pattern = "RNA", invert = T)] ),
                          selected = "None")
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
                                       cols = c("lightgrey","firebrick1"),
                                       min.cutoff = 0,
                                       reduction = input$red,
                                       label = TRUE,
                                       combine = FALSE))
                })
            })
        }
    })
    
    output$out_feat <- renderUI({
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
                cat = input$split
                sp = FALSE
                ii <- i
                output[[paste0('plot_ridge',ii)]] <- renderPlot({
                    if ( input$split.plot == "Yes" ) { 
                        sp = TRUE 
                    } else { sp = FALSE }
                    if ( input$split == "None" ) {
                        cat = NULL
                        sp = FALSE
                    }
                    if ( input$type == "RidgePlot" ) {
                        return(RidgePlot(datasetInput(),
                                         features = input$genes[[ii]],
                                         combine = FALSE))
                    }
                    if ( input$type == "VlnPlot" ) {
                        return(VlnPlot(datasetInput(), 
                                       pt.size = 0,
                                       split.by = cat,
                                       split.plot = sp,
                                       features = input$genes[[ii]],
                                       combine = FALSE))
                    }
                    if ( input$type == "DotPlot" ) {
                        return(DotPlot(datasetInput(), 
                                       features = input$genes[[ii]],
                                       cols = "RdYlBu"))
                    }
                    if (input$type == "HeatMap") {
                        return(DoHeatmap(datasetInput(),
                                         features = input$genes[[ii]],
                                         draw.lines = F,
                                         size = 4))
                    }
                })
            })
        }
    })
    
    output$out_multi <- renderUI({
        if (length(input$genes) <= 1){return(NULL)}
        out <- plotOutput(outputId = "plot_multi")
        return(out)
    })
    
    observe({
        output[["plot_multi"]] <- renderPlot({
            if (input$type == "RidgePlot") {
                return(RidgePlot(datasetInput(),
                                 features = input$genes,
                                 stack = T))
            }
            if (input$type == "VlnPlot") {
                cat = input$split
                sp = FALSE
                if ( input$split.plot == "Yes" ) { 
                    sp = TRUE 
                } else { sp = FALSE }
                if ( input$split == "None" ) {
                    cat = NULL
                    sp = FALSE
                }
                return(VlnPlot(datasetInput(),
                               features = input$genes,
                               stack = T,
                               flip = T,
                               split.by = cat,
                               split.plot = sp))
            }
            if (input$type == "DotPlot") {
                return(DotPlot(datasetInput(),
                               features = input$genes,
                               cols = "RdYlBu"))
            }
            if (input$type == "HeatMap") {
                return(DoHeatmap(datasetInput(),
                                 features = input$genes,
                                 size = 3.5))
            }
        })
    })
    
    output$out_multi2 <- renderUI({
        if (length(input$genes) <= 1) {return(NULL)}
        if (input$type == "VlnPlot") {
            out <- plotOutput(outputId = "plot_multi2")
        } else { out <- NULL }
        return(out)
    })
    
    observe({
        output[["plot_multi2"]] <- renderPlot({
            cat = input$split
            sp = FALSE
            if ( input$split.plot == "Yes" ) { 
                sp = TRUE 
            } else { sp = FALSE }
            if ( input$split == "None" ) {
                cat = NULL
                sp = FALSE
            }
            return(VlnPlot(datasetInput(),
                           features = input$genes,
                           split.by = cat,
                           split.plot = sp))
        })
    })
    
    output$out_table <- renderTable({ return(Avg())} , rownames = T)
    output$export_table <- downloadHandler(
        filename = function() { paste("AverageExpression.csv") },
        content = function(file) { write.csv(Avg(), file, row.names = TRUE, quote = FALSE, col.names = NULL) }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
