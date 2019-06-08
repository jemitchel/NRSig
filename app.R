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
                   label = "Compute"),
      # downloadButton("downloadfRMA",label = "Download fRMA Preprocessed Data"),
      # downloadButton("downloadDiffex",label = "Download Differential Expression Data")
      uiOutput("download1"),
      uiOutput("download2"),
      uiOutput("download3")
    ),

    mainPanel(
      DT::dataTableOutput("resTable"),
      plotOutput("bplots")
      # downloadButton("downloadPlots",label = "Download Plot")
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
  
  # produces the results dataframe with column of buttons
  tmptbl <- reactive({data.frame(res()[[2]],
                      Actions = shinyInput(actionButton, 15, 'button_', label = "See Targets", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' )
  )})
  
  # renders the results dataframe with column of buttons
  output$resTable <- DT::renderDataTable({
    tmptbl()
  },server = FALSE, escape = FALSE, selection = 'none')

  # chooses the correct plot depending which button is pushed
  chosenPlt <- eventReactive(input$select_button, {
    selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    selectedNR <- rownames(tmptbl())[selectedRow]
    return(res()[[1]][[selectedNR]])
  })

  # shows the selected table
  observe({
    output$bplots <- renderPlot({
      chosenPlt()[[1]]
    },height = (5000/74)*chosenPlt()[[2]],width = 500)
  })
  
  
  output$download1 <- renderUI({
    if(!is.null(res())) {
      downloadButton('downloadfRMA', label = "Download fRMA Preprocessed Data")
    }
  })
  
  output$download2 <- renderUI({
    if(!is.null(res())) {
      downloadButton('downloadDiffex', label = "Download Differential Expression Data")
    }
  })
    
  output$download3 <- renderUI({
    if(!is.null(chosenPlt())) {
      downloadButton('downloadPlot', label = "Download Plot")
    }
  })
  
  output$downloadfRMA <- downloadHandler(
    filename = function() { paste('input_preprocessed', '.csv', sep='') },
    content = function(file) {
      write.csv("INSERT TABLE HERE", file) #need to make the table first
    }
  )
  
  output$downloadDiffex <- downloadHandler(
    filename = function() { paste('diff_exprs', '.csv', sep='') },
    content = function(file) {
      write.csv("INSERT TABLE HERE", file) #need to make the table first
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(chosenPlt()[[1]]$labels$title, '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = chosenPlt()[[1]], device = "png",
             height = (60/74)*chosenPlt()[[2]], width = 5)
    }
  )
  
}

shinyApp(ui, server)






