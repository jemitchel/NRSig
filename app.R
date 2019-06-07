library(shiny)
library(DT)

ui <- fluidPage(

  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = "files",
                label = "Upload .CEL files for one sample or multiple replicates (replicates will be averaged)",
                multiple = TRUE),
      verbatimTextOutput("fileName"),
      actionButton(inputId = "compute",
                   label = "Compute")
    ),
    mainPanel(
      DT::dataTableOutput("resTable"),
      plotOutput("bplots")
    )
  )
  
  # fluidRow(
  #   column(3,
  #          fileInput(inputId = "files",
  #                    label = "Upload .CEL files for one sample or multiple replicates (replicates will be averaged)",
  #                    multiple = TRUE),
  #          verbatimTextOutput("fileName"),
  #          actionButton(inputId = "compute",
  #                       label = "Compute")),
  #   column(3,
  #          DT::dataTableOutput("resTable")),
  #   column(6,
  #          plotOutput("bplots"))
  # )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)

  output$fileName <- renderText({
    tmp <- "Uploaded Files:\n"
    if (!is.null(input$files$datapath)) {
      for (i in basename(input$files$name)) {
        tmp <- paste(tmp,i,"\n")
      }
      return(tmp)
    } else {
      return(NULL)
    }
  })


  res <- eventReactive(input$compute, {
    progress <- Progress$new(session)
    progress$set(message = 'Preprocessing',
                 detail = 'This may take a few minutes...')
    source("C:/Users/jonat/Documents/R/NRSig-app/preprocess.R")
    samples_matrix <- pre_proc(input$files$datapath)
    progress$set(message = 'Preprocessing Completed', detail = "")
    on.exit(progress$close())

    withProgress(message = 'Computing Enrichment of NR-Target Gene Sets...', value = 0,{
      source("C:/Users/jonat/Documents/R/NRSig-app/compute_enriched_NRs.R")
      results <- CalcEnrich(samples_matrix)
      return(results)
    })

  })

  # makes inputs for buttons on output results table
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }

  # generates the results dataframe with column of buttons
  output$resTable <- DT::renderDataTable({data.frame(res()[[2]],
    Actions = shinyInput(actionButton, 15, 'button_', label = "See Targets", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' )
  )},server = FALSE, escape = FALSE, selection = 'none')

  # chooses the correct plot depending which button is pushed
  chosenPlt <- eventReactive(input$select_button, {
    selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    return(res()[[1]][[selectedRow]])
  })

  # shows the selected table
  observe({
    output$bplots <- renderPlot({
      chosenPlt()[[1]]
    },height = (5000/74)*chosenPlt()[[2]],width = 500)
  })
  

}

shinyApp(ui, server)










# library(shiny)
# library(DT)
# 
# shinyApp(
#   ui <- fluidPage(
#     DT::dataTableOutput("data"),
#     DT::dataTableOutput("tbl")
#   ),
# 
#   server <- function(input, output) {
# 
#     shinyInput <- function(FUN, len, id, ...) {
#       inputs <- character(len)
#       for (i in seq_len(len)) {
#         inputs[i] <- as.character(FUN(paste0(id, i), ...))
#       }
#       inputs
#     }
# 
#     # generates the dataframe
#     output$data <- DT::renderDataTable({data.frame(Name = c('Dilbert', 'Alice', 'Wally', 'Ashok', 'Dogbert'),
#       Actions = shinyInput(actionButton, 5, 'button_', label = "See Targets", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' )
#     )},server = FALSE, escape = FALSE, selection = 'none')
# 
#     # # generates the dataframe
#     # tmp <- reactive({data.frame(Name = c('Dilbert', 'Alice', 'Wally', 'Ashok', 'Dogbert'),
#     #                                                Actions = shinyInput(actionButton, 5, 'button_', label = "See Targets", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' )
#     # )},server = FALSE, escape = FALSE, selection = 'none')
#     # 
#     # output$data <- DT::renderDataTable({})
#     #   
#     
#     # df <- reactiveValues(data = data.frame(
#     #
#     #   Name = c('Dilbert', 'Alice', 'Wally', 'Ashok', 'Dogbert'),
#     #   Motivation = c(62, 73, 3, 99, 52),
#     #   Actions = shinyInput(actionButton, 5, 'button_', label = "Fire", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),
#     #   stringsAsFactors = FALSE,
#     #   row.names = 1:5
#     # ))
#     #
#     #
#     # output$data <- DT::renderDataTable(
#     #   df$data, server = FALSE, escape = FALSE, selection = 'none'
#     # )
#     #
#     oot <- eventReactive(input$select_button, {
#       selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
#       return(clst[[selectedRow]])
#     })
# 
# 
#     output$tbl <- DT::renderDataTable({
#       oot()
#     })
# 
#   }
# )





